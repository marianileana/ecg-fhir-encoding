"""Step 1: QRS detection on the gateway.

The Pan-Tompkins detector (NeuroKit2 implementation) is applied as a gateway would apply
it: in 60-s blocks, with 2 s of signal before and 1 s after each block. The detected R
peaks are compared with the reference beat annotations (tolerance 150 ms, ANSI/AAMI EC57).

Output: results/detection.json (Se, PPV, TP/FN/FP per record and in total, processor
time) and results/peaks_<record>.npy (R peaks, used by the next scripts).
"""
import json
import os
import time

import neurokit2 as nk
import numpy as np
from wfdb import processing

from common import FS, OUT, load, records

BLOCK, BACK, AHEAD = 60 * FS, 2 * FS, 1 * FS   # block length, look-back and look-ahead (samples)
TOL = int(round(0.150 * FS))                    # 150 ms = 54 samples


def detect_stream(sig_mv):
    """R-peak positions (samples) for a whole record, processed block by block."""
    n = len(sig_mv)
    peaks = []
    for start in range(0, n, BLOCK):
        lo, hi = max(0, start - BACK), min(n, start + BLOCK + AHEAD)
        seg = sig_mv[lo:hi]
        clean = nk.ecg_clean(seg, sampling_rate=FS, method="pantompkins1985")
        _, info = nk.ecg_peaks(clean, sampling_rate=FS, method="pantompkins1985")
        p = np.asarray(info["ECG_R_Peaks"], dtype=np.int64) + lo
        p = p[(p >= start) & (p < min(n, start + BLOCK))]   # keep only peaks inside the block
        for q in p:
            # a detection within 150 ms of the previous one is a duplicate at a block border
            if not peaks or q - peaks[-1] > TOL:
                peaks.append(int(q))
    return np.array(peaks, dtype=np.int64)


def main():
    per = {}
    tot_t = 0.0
    tot_s = 0
    for rec in records():
        sig, lead, gain, base, ref, sym = load(rec)
        mv = (sig - base) / gain
        t0 = time.process_time()
        det = detect_stream(mv)
        dt = time.process_time() - t0
        # wfdb counts a match only when |difference| < window, so TOL + 1 gives the inclusive 150 ms window
        cmp = processing.compare_annotations(ref, det, TOL + 1)
        tp, fn, fp = int(cmp.tp), int(cmp.fn), int(cmp.fp)
        per[rec] = dict(lead=lead, n_ref=int(len(ref)), n_det=int(len(det)), tp=tp, fn=fn, fp=fp,
                        se=tp / (tp + fn), ppv=tp / (tp + fp), cpu_s=dt, n_samples=int(len(sig)))
        np.save(os.path.join(OUT, f"peaks_{rec}.npy"), det)
        tot_t += dt
        tot_s += len(sig)
        print(f"{rec} {lead:5s} ref={len(ref):5d} det={len(det):5d} Se={per[rec]['se']:.4f} "
              f"PPV={per[rec]['ppv']:.4f} cpu={dt:.2f}s", flush=True)

    TP = sum(v["tp"] for v in per.values())
    FN = sum(v["fn"] for v in per.values())
    FP = sum(v["fp"] for v in per.values())
    hours = tot_s / FS / 3600
    summary = dict(records=len(per), beats_ref=TP + FN, tp=TP, fn=FN, fp=FP,
                   se=TP / (TP + FN), ppv=TP / (TP + FP),
                   f1=2 * TP / (2 * TP + FP + FN),
                   hours=hours, cpu_s_total=tot_t, cpu_s_per_hour=tot_t / hours,
                   tolerance_ms=150, block_s=60, lookback_s=2, lookahead_s=1,
                   detector="Pan-Tompkins (NeuroKit2 0.2.13 implementation, method=pantompkins1985)")
    print(json.dumps(summary, indent=1))
    with open(os.path.join(OUT, "detection.json"), "w") as f:
        json.dump(dict(summary=summary, per_record=per), f, indent=1)


if __name__ == "__main__":
    main()
