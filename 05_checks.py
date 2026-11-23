"""Step 5: checks behind statements in Sections 3 and 4 of the article.

- negative control of the FHIR validation: invalid resources must be rejected (Section 4.2)
- lead used in each record (Section 4.1)
- integrity of the database against its SHA256SUMS.txt (Section 4.1)
- SampledData parameters and converter gain and baseline of every record (Section 3.2)
- record 100 compared with the copy distributed with the wfdb software (Section 4.1);
  set WFDB_SAMPLE to the wfdb-python "sample-data" folder (default: wfdbrepo/sample-data)

Output: results/checks.json.
"""
import hashlib
import json
import os

from fhir.resources.R4B.observation import Observation

from common import DB, HERE, OUT, ecg_observation, hr_observation, load, records


def main():
    out = {}

    # negative control of the validation
    sig, lead, gain, base, _, _ = load("100")
    good = ecg_observation("100", lead, sig[:360], 0, gain, base)
    cases = {
        "valid ECG Observation": good,
        "valid HR Observation": hr_observation("100", 0, 10, 72.0),
        "missing status": {k: v for k, v in good.items() if k != "status"},
        "non-numeric period": {**good, "valueSampledData": {**good["valueSampledData"], "period": "x"}},
        "unknown element": {**good, "foo": 1},
    }
    res = {}
    for name, obj in cases.items():
        try:
            Observation.model_validate(obj)
            res[name] = "accepted"
        except Exception:
            res[name] = "rejected"
    out["validation_control"] = res

    # lead used in each record
    leads = {r: load(r)[1] for r in records()}
    out["lead_counts"] = {k: sum(1 for v in leads.values() if v == k) for k in sorted(set(leads.values()))}
    out["non_MLII_records"] = sorted(r for r, v in leads.items() if v != "MLII")

    # integrity of the database files
    ok = bad = 0
    for line in open(os.path.join(DB, "SHA256SUMS.txt")):
        h, name = line.split()
        with open(os.path.join(DB, name), "rb") as f:
            if hashlib.sha256(f.read()).hexdigest() == h:
                ok += 1
            else:
                bad += 1
    out["sha256"] = {"ok": ok, "failed": bad}

    # SampledData parameters (record 100)
    sd = good["valueSampledData"]
    out["sampleddata"] = {"gain": gain, "baseline": base, "origin_mV": sd["origin"]["value"],
                          "factor": sd["factor"], "period_ms": sd["period"]}

    # gain and baseline of every record
    gb = {}
    for r in records():
        _, _, g, b0, _, _ = load(r)
        gb.setdefault(f"{g:g}/{b0}", []).append(r)
    out["gain_baseline"] = {k: len(v) for k, v in gb.items()}

    # record 100 against the copy distributed with wfdb-python
    sample = os.environ.get("WFDB_SAMPLE", os.path.join(HERE, "wfdbrepo", "sample-data"))
    cmp100 = {}
    for ext in ("dat", "atr"):
        a, b = os.path.join(DB, f"100.{ext}"), os.path.join(sample, f"100.{ext}")
        if os.path.exists(b):
            with open(a, "rb") as fa, open(b, "rb") as fb:
                cmp100[ext] = fa.read() == fb.read()
        else:
            cmp100[ext] = "copy not found"
    out["record100_vs_wfdb_sample"] = cmp100

    with open(os.path.join(OUT, "checks.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
