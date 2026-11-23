"""Step 3: additional values reported in Sections 5.3 and 5.4 of the article.

- S5 payload per record and share of 10-s windows with beats outside class N (Fig. 4)
- records without ectopic beats and their S5 payload
- false-alarm windows that contain no beat annotation (e.g. ventricular flutter in record 207)
- QRS detection on whole records, to measure the effect of the 60-s blocks

Input: results/peaks_<record>.npy from 01_detect.py. Output: results/extra.json.
"""
import importlib
import json
import os
import warnings

import neurokit2 as nk
import numpy as np
from wfdb import processing

from common import AAMI, FS, OUT, PACED, load, records

warnings.filterwarnings("ignore")
enc = importlib.import_module("02_encode_send")   # reuse the gateway logic of step 2
W = enc.W


def main():
    per = {}
    tp = fn = fp = 0
    ectopic_w = all_w = 0
    false_empty = {}
    for rec in records():
        sig, lead, gain, base, ref, sym = load(rec)
        n = len(sig)
        hours = n / FS / 3600
        peaks = np.load(os.path.join(OUT, f"peaks_{rec}.npy"))
        hr = enc.hr_windows(peaks, n)
        s5, trig = enc.build_s5(rec, sig, lead, gain, base, peaks, hr, enc.THETA_MAIN)
        cls = np.array([AAMI.get(s, "Q") for s in sym])
        nwin = (n - 1) // W + 1
        ect_w = set((ref[cls != "N"] // W).tolist())
        per[rec] = dict(paced=rec in PACED, s5_kB_per_hour=sum(len(b) for _, b, _ in s5) / hours / 1000,
                        trig_pct=100 * len(trig) / nwin, ectopic_window_pct=100 * len(ect_w) / nwin,
                        n_veb=int(np.sum(cls == "VEB")), n_sveb=int(np.sum(cls == "SVEB")))
        if rec not in PACED:
            ectopic_w += len(ect_w)
            all_w += nwin
            beat_w = set((ref // W).tolist())
            empty = sorted(int(w) for w in trig - ect_w if w not in beat_w)
            if empty:   # sent windows without any beat annotation
                false_empty[rec] = empty

        # whole-record detection (no 60-s blocks); 55 = inclusive 150 ms window (54 samples)
        mv = (sig - base) / gain
        clean = nk.ecg_clean(mv, sampling_rate=FS, method="pantompkins1985")
        _, info = nk.ecg_peaks(clean, sampling_rate=FS, method="pantompkins1985")
        c = processing.compare_annotations(ref, np.asarray(info["ECG_R_Peaks"], dtype=int), 55)
        tp += c.tp
        fn += c.fn
        fp += c.fp

    v = np.array([p["s5_kB_per_hour"] for p in per.values()])
    quiet = [r for r, p in per.items() if p["n_veb"] + p["n_sveb"] == 0 and not p["paced"]]
    out = dict(
        per_record=per,
        s5_kB_per_hour=dict(min=float(v.min()), median=float(np.median(v)), max=float(v.max()),
                            argmin=min(per, key=lambda r: per[r]["s5_kB_per_hour"]),
                            argmax=max(per, key=lambda r: per[r]["s5_kB_per_hour"])),
        records_without_ectopy=quiet,
        s5_quiet_median_kB_per_hour=float(np.median([per[r]["s5_kB_per_hour"] for r in quiet])) if quiet else None,
        oracle_ectopic_window_pct=100 * ectopic_w / all_w,
        false_windows_without_beats=dict(count=sum(len(x) for x in false_empty.values()), windows=false_empty),
        whole_record_detection=dict(tp=int(tp), fn=int(fn), fp=int(fp), se=tp / (tp + fn), ppv=tp / (tp + fp)),
    )
    with open(os.path.join(OUT, "extra.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps({k: x for k, x in out.items() if k != "per_record"}, indent=1))


if __name__ == "__main__":
    main()
