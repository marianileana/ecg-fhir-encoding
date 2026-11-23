"""Step 2: FHIR encodings, transmission to a FHIR server and evaluation.

Encodings (one single-lead ECG per record):
    S1   raw ECG, one Observation per 1 s, one request per Observation (JSON)
    S1b  raw ECG, one Observation per 10 s, one request per Observation (JSON);
         S1 -> S1b changes only the length of the Observation,
         S1b -> S2 changes only the grouping in Bundles
    S2   raw ECG, transaction Bundle of six 10-s Observations per 60 s (JSON)
    S3   as S2, with gzip
    S4   heart rate only, Bundle of 10-s HR Observations per 60 s, gzip
    S5   hybrid mode: S4 plus the 10-s ECG segments that contain a premature beat, gzip

The FHIR server runs in a separate process on the same computer and is reached over
HTTP/1.1 (loopback); a single client sends the requests one after another. The server
decompresses, parses and validates every request but does not store the resources.

Also computed: heart rate error against the annotations, share of ectopic beats kept by
the hybrid mode for several thresholds, and a check that the ECG mapping is lossless.

Input: results/peaks_<record>.npy and results/detection.json from 01_detect.py.
Output: results/encode_send.json.
"""
import gzip
import http.client
import json
import multiprocessing as mp
import os
import socket
import time
import zlib
from http.server import BaseHTTPRequestHandler, HTTPServer

import numpy as np

from common import (FS, OUT, PACED, AAMI, bundle, dumps, ecg_observation, flag_interpretation,
                    hr_observation, load, records)

PORT = 18080
W = 10 * FS                 # 10-s window (samples)
MIN = 60 * FS               # one Bundle per 60 s (samples)
THETAS = [0.70, 0.75, 0.80, 0.85, 0.90, 0.95]   # prematurity thresholds for the sweep (Fig. 3)
THETA_MAIN = 0.85           # threshold fixed before the experiment
MED_N = 8                   # median of the previous 8 RR intervals
GZIP_LEVEL = 6
STRATEGIES = ["S1", "S1b", "S2", "S3", "S4", "S5"]


# ======================================================================= FHIR server
def serve(port, ready):
    """Minimal FHIR endpoint: decompress, parse and validate (FHIR R4B models), no storage.
    Replies with the validation result and the server processing time."""
    from fhir.resources.R4B.bundle import Bundle
    from fhir.resources.R4B.observation import Observation

    class Handler(BaseHTTPRequestHandler):
        protocol_version = "HTTP/1.1"       # persistent connections
        disable_nagle_algorithm = True      # TCP_NODELAY on the server side

        def log_message(self, *args):
            pass

        def do_POST(self):
            n = int(self.headers["Content-Length"])
            body = self.rfile.read(n)
            t0 = time.perf_counter()
            if self.headers.get("Content-Encoding") == "gzip":
                body = gzip.decompress(body)
            obj = json.loads(body)
            ok, nres, err = True, 0, ""
            try:
                if obj.get("resourceType") == "Bundle":
                    b = Bundle.model_validate(obj)
                    nres = len(b.entry or [])
                else:
                    Observation.model_validate(obj)
                    nres = 1
            except Exception as e:
                ok, err = False, str(e)[:200]
            dt = (time.perf_counter() - t0) * 1000
            out = json.dumps({"ok": ok, "n": nres, "ms": dt, "err": err}).encode()
            self.send_response(200 if ok else 422)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(out)))
            self.end_headers()
            self.wfile.write(out)

    srv = HTTPServer(("127.0.0.1", port), Handler)
    ready.set()
    srv.serve_forever()


class Client:
    """HTTP/1.1 client with one persistent connection. Counts every byte written to the
    socket (request line, headers and body) and measures the round-trip time."""

    def __init__(self, port):
        self.c = http.client.HTTPConnection("127.0.0.1", port)
        self.c.connect()
        self.c.sock.setsockopt(socket.IPPROTO_TCP, socket.TCP_NODELAY, 1)
        self.sent = 0
        orig = self.c.send

        def counting_send(data):
            self.sent += len(data)
            return orig(data)

        self.c.send = counting_send

    def post(self, path, body, gz):
        """Send one request; return (bytes sent, round-trip time in ms, server reply)."""
        h = {"Content-Type": "application/fhir+json"}
        if gz:
            h["Content-Encoding"] = "gzip"
        before = self.sent
        t0 = time.perf_counter()
        self.c.request("POST", path, body=body, headers=h)
        r = self.c.getresponse()
        data = json.loads(r.read())
        lat = (time.perf_counter() - t0) * 1000
        return self.sent - before, lat, data


# ===================================================================== gateway logic
def hr_windows(peaks, n):
    """Heart rate (bpm) per 10-s window: 60 / median RR (s) of the RR intervals that end in
    the window; at least two intervals are required. Returns {window index: bpm}."""
    out = {}
    if len(peaks) < 2:
        return out
    rr = np.diff(peaks)
    w = peaks[1:] // W          # window of the later beat of each interval
    for k in np.unique(w):
        r = rr[w == k]
        if len(r) >= 2:
            out[int(k)] = 60.0 * FS / float(np.median(r))
    return out


def triggers(peaks, theta):
    """Premature-beat trigger (eq. (2)): beat i is premature when
    RR_i < theta * median(RR_{i-8}, ..., RR_{i-1}).
    Returns the set of 10-s windows that contain a premature beat."""
    rr = np.diff(peaks)
    trig = set()
    for i in range(MED_N, len(rr)):
        med = np.median(rr[i - MED_N:i])
        if rr[i] < theta * med:
            trig.add(int(peaks[i + 1] // W))
    return trig


def build(rec, sig, lead, gain, base, peaks, theta):
    """Build the messages of all strategies for one record.

    Returns ({strategy: [(path, body, gzip?)]}, {strategy: gateway processor time (s)},
    trigger windows, heart rate per window). The QRS detection time is added later."""
    n = len(sig)
    msgs, cpu = {}, {}

    # S1: one 1-s Observation per request
    t = time.process_time()
    m = []
    for s in range(0, n, FS):
        o = ecg_observation(rec, lead, sig[s:s + FS], s / FS, gain, base)
        m.append(("/fhir/Observation", dumps(o).encode(), False))
    cpu["S1"] = time.process_time() - t
    msgs["S1"] = m

    # S1b: the 10-s Observations of S2, one per request
    t = time.process_time()
    m = []
    for s in range(0, n, W):
        o = ecg_observation(rec, lead, sig[s:min(s + W, n)], s / FS, gain, base)
        m.append(("/fhir/Observation", dumps(o).encode(), False))
    cpu["S1b"] = time.process_time() - t
    msgs["S1b"] = m

    # S2: six 10-s Observations per transaction Bundle
    t = time.process_time()
    raw_bundles = []
    for b0 in range(0, n, MIN):
        obs = [ecg_observation(rec, lead, sig[s:min(s + W, n)], s / FS, gain, base)
               for s in range(b0, min(b0 + MIN, n), W)]
        raw_bundles.append(dumps(bundle(obs, f"{rec}/raw/{b0}")).encode())
    cpu["S2"] = time.process_time() - t
    msgs["S2"] = [("/fhir", b, False) for b in raw_bundles]

    # S3: the S2 Bundles compressed with gzip (mtime=0 keeps the output identical between runs)
    t = time.process_time()
    msgs["S3"] = [("/fhir", gzip.compress(b, GZIP_LEVEL, mtime=0), True) for b in raw_bundles]
    cpu["S3"] = cpu["S2"] + (time.process_time() - t)

    # S4: heart rate only, one Bundle per minute, gzip
    t = time.process_time()
    hr = hr_windows(peaks, n)
    cpu_hr = time.process_time() - t
    s4 = []
    t = time.process_time()
    for b0 in range(0, n, MIN):
        ws = range(b0 // W, (min(b0 + MIN, n) - 1) // W + 1)
        hro = [hr_observation(rec, k * 10, 10, hr[k]) for k in ws if k in hr]
        if hro:
            s4.append(("/fhir", gzip.compress(dumps(bundle(hro, f"{rec}/hr/{b0}")).encode(),
                                              GZIP_LEVEL, mtime=0), True))
    cpu["S4"] = cpu_hr + (time.process_time() - t)
    msgs["S4"] = s4

    # S5: hybrid mode (includes the trigger)
    t = time.process_time()
    s5, trig = build_s5(rec, sig, lead, gain, base, peaks, hr, theta)
    cpu["S5"] = cpu_hr + (time.process_time() - t)
    msgs["S5"] = s5
    return msgs, cpu, trig, hr


def build_s5(rec, sig, lead, gain, base, peaks, hr, theta):
    """Hybrid mode: per minute, the HR Observations plus the 10-s ECG segments of the
    windows marked by the trigger (with the Abnormal interpretation code), gzip."""
    n = len(sig)
    s5 = []
    trig = triggers(peaks, theta)
    for b0 in range(0, n, MIN):
        ws = range(b0 // W, (min(b0 + MIN, n) - 1) // W + 1)
        res = [hr_observation(rec, k * 10, 10, hr[k]) for k in ws if k in hr]
        res += [flag_interpretation(ecg_observation(rec, lead, sig[k * W:min((k + 1) * W, n)], k * 10,
                                                    gain, base)) for k in ws if k in trig]
        if res:
            s5.append(("/fhir", gzip.compress(dumps(bundle(res, f"{rec}/hyb/{b0}")).encode(),
                                              GZIP_LEVEL, mtime=0), True))
    return s5, trig


def check_lossless(msgs, sig, gain, base):
    """Decode every ECG resource of the messages and compare the values in mV with the
    database values. Returns (largest absolute error, number of resources checked)."""
    ref = (sig - base) / gain
    maxerr, nchk = 0.0, 0
    for _, body, gz in msgs:
        obj = json.loads(gzip.decompress(body) if gz else body)
        res = [e["resource"] for e in obj["entry"]] if obj["resourceType"] == "Bundle" else [obj]
        for r in res:
            sd = r.get("valueSampledData")
            if not sd:
                continue
            d = np.array(sd["data"].split(), dtype=np.float64)
            mv = sd["origin"]["value"] + sd["factor"] * d
            # sample index of the first value, from the timestamp of the resource
            s = int(round((np.datetime64(r["effectiveDateTime"][:-1]) -
                           np.datetime64("2026-01-01T00:00:00")) / np.timedelta64(1, "ms") * FS / 1000))
            maxerr = max(maxerr, float(np.max(np.abs(mv - ref[s:s + len(d)]))))
            nchk += 1
    return maxerr, nchk


# ============================================================================= main
def main():
    # start the FHIR server in its own process
    ready = mp.Event()
    p = mp.Process(target=serve, args=(PORT, ready), daemon=True)
    p.start()
    ready.wait()
    cli = Client(PORT)

    det = json.load(open(os.path.join(OUT, "detection.json")))["per_record"]
    agg = {k: dict(msgs=0, body=0, wire=0, resources=0, valid_msgs=0, invalid_msgs=0,
                   cpu_s=0.0, lat=[], val_ms=[]) for k in STRATEGIES}
    hr_err = []
    hr_windows_total, hr_windows_edge = 0, 0
    sweep = {th: dict(trig_windows=0, windows=0, samples=0, veb=0, veb_cap=0, sveb=0, sveb_cap=0,
                      false_windows=0, body=0) for th in THETAS}
    lossless = dict(maxerr=0.0, resources=0)
    total_samples = 0

    for rec in records():
        sig, lead, gain, base, ref, sym = load(rec)
        n = len(sig)
        total_samples += n
        peaks = np.load(os.path.join(OUT, f"peaks_{rec}.npy"))
        msgs, cpu, trig, hr = build(rec, sig, lead, gain, base, peaks, THETA_MAIN)
        # S4 and S5 also need the QRS detection on the gateway
        cpu["S4"] += det[rec]["cpu_s"]
        cpu["S5"] += det[rec]["cpu_s"]

        # send everything and collect bytes, times and validation results
        for k in STRATEGIES:
            a = agg[k]
            for path, body, gz in msgs[k]:
                wire, lat, resp = cli.post(path, body, gz)
                a["msgs"] += 1
                a["body"] += len(body)
                a["wire"] += wire
                a["resources"] += resp["n"]
                a["valid_msgs" if resp["ok"] else "invalid_msgs"] += 1
                if not resp["ok"]:
                    print("INVALID", rec, k, resp["err"])
                a["lat"].append(lat)
                a["val_ms"].append(resp["ms"])
            a["cpu_s"] += cpu[k]

        # lossless check on the ECG resources of S1 and S3
        for k in ("S1", "S3"):
            e, c = check_lossless(msgs[k], sig, gain, base)
            lossless["maxerr"] = max(lossless["maxerr"], e)
            lossless["resources"] += c

        # heart rate error: gateway HR against the HR computed in the same way from the annotations
        hr_ref = hr_windows(ref, n)
        hr_windows_total += len(hr_ref)
        hr_windows_edge += len(hr)
        for w_, v in hr_ref.items():
            if w_ in hr:
                hr_err.append(hr[w_] - v)

        # ectopic beats kept by the hybrid mode, non-paced records only (AAMI practice)
        if rec not in PACED:
            cls = np.array([AAMI.get(s, "Q") for s in sym])
            nwin = (n - 1) // W + 1
            ect_w = set((ref[cls != "N"] // W).tolist())   # windows with a beat outside class N
            for th in THETAS:
                tg = trig if th == THETA_MAIN else triggers(peaks, th)
                sw = sweep[th]
                sw["windows"] += nwin
                sw["samples"] += n
                sw["trig_windows"] += len(tg)
                sw["false_windows"] += len(tg - ect_w)     # sent windows without non-N beats
                for c_ in ("VEB", "SVEB"):
                    idx = ref[cls == c_] // W
                    sw[c_.lower()] += len(idx)
                    sw[c_.lower() + "_cap"] += int(sum(1 for w_ in idx if int(w_) in tg))
                if th == THETA_MAIN:
                    sw["body"] += sum(len(b) for _, b, _ in msgs["S5"])
                else:
                    m2, _ = build_s5(rec, sig, lead, gain, base, peaks, hr, th)
                    sw["body"] += sum(len(b) for _, b, _ in m2)
        print(f"{rec} done  S1={len(msgs['S1'])} S5={len(msgs['S5'])} trig={len(trig)}", flush=True)

    # ------------------------------------------------------------- summary per hour of ECG
    hours = total_samples / FS / 3600
    summary = {}
    for k in STRATEGIES:
        a = agg[k]
        lat = np.array(a["lat"])
        v = np.array(a["val_ms"])
        summary[k] = dict(
            msgs=a["msgs"], msgs_per_hour=a["msgs"] / hours,
            resources=a["resources"], valid_msgs=a["valid_msgs"], invalid_msgs=a["invalid_msgs"],
            body_bytes=a["body"], wire_bytes=a["wire"],
            body_kB_per_hour=a["body"] / hours / 1000, wire_kB_per_hour=a["wire"] / hours / 1000,
            wire_kbps=a["wire"] * 8 / (hours * 3600) / 1000,
            edge_cpu_s_per_hour=a["cpu_s"] / hours,
            lat_ms_median=float(np.median(lat)), lat_ms_p95=float(np.percentile(lat, 95)),
            val_ms_median=float(np.median(v)), val_ms_p95=float(np.percentile(v, 95)),
            ingest_ms_per_hour=float(lat.sum() / hours))
    base_wire = summary["S1"]["wire_bytes"]
    for k in STRATEGIES:
        summary[k]["wire_reduction_vs_S1_pct"] = 100 * (1 - summary[k]["wire_bytes"] / base_wire)

    raw_bytes_11bit = total_samples * 11 / 8   # samples packed at 11 bits, the reference size
    e = np.array(hr_err)
    hrsum = dict(windows_ref=hr_windows_total, windows_edge=hr_windows_edge, pairs=int(len(e)),
                 mae_bpm=float(np.mean(np.abs(e))), median_abs_bpm=float(np.median(np.abs(e))),
                 p95_abs_bpm=float(np.percentile(np.abs(e), 95)),
                 within_5bpm_pct=float(100 * np.mean(np.abs(e) <= 5)), bias_bpm=float(np.mean(e)))
    sw_out = {}
    for th, s in sweep.items():
        sw_hours = s["samples"] / FS / 3600
        sw_out[f"{th:.2f}"] = dict(
            **s, veb_se=s["veb_cap"] / s["veb"], sveb_se=s["sveb_cap"] / s["sveb"],
            raw_time_pct=100 * s["trig_windows"] / s["windows"],
            false_windows_per_hour=s["false_windows"] / sw_hours,
            body_kB_per_hour=s["body"] / sw_hours / 1000)
    out = dict(hours=hours, records=len(records()), samples=total_samples,
               raw_11bit_kB_per_hour=raw_bytes_11bit / hours / 1000,
               strategies=summary, hr=hrsum, sweep=sw_out, theta_main=THETA_MAIN,
               lossless=lossless, fhir="R4 resources validated with fhir.resources R4B models",
               zlib=zlib.ZLIB_RUNTIME_VERSION)
    with open(os.path.join(OUT, "encode_send.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps({k: {kk: round(vv, 3) if isinstance(vv, float) else vv for kk, vv in v.items()}
                      for k, v in summary.items()}, indent=1))
    print(json.dumps(hrsum, indent=1))
    print(json.dumps(sw_out, indent=1))
    print("lossless", lossless)
    p.terminate()


if __name__ == "__main__":
    main()
