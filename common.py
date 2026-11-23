"""Shared helpers used by all scripts: paths, data loading, AAMI beat classes and the
builders of the HL7 FHIR resources sent by the gateway.

Data: MIT-BIH Arrhythmia Database v1.0.0 (PhysioNet), 48 half-hour records sampled at
360 Hz with 11-bit resolution; converter gain 200 adu/mV and baseline 1024 in all records.
"""
import json
import os
import uuid
from datetime import datetime, timedelta, timezone

import numpy as np
import wfdb

# ------------------------------------------------------------------ paths and constants
HERE = os.path.dirname(os.path.abspath(__file__))
DB = os.environ.get("MITDB", os.path.join(HERE, "data", "mit-bih-arrhythmia-database-1.0.0"))
OUT = os.path.join(HERE, "results")
os.makedirs(OUT, exist_ok=True)

FS = 360                                   # sampling frequency (Hz)
PACED = {"102", "104", "107", "217"}       # records with paced beats, excluded from the ectopy analysis

# MIT-BIH beat annotation symbols and their AAMI EC57 classes
# (N: normal, SVEB: supraventricular ectopic, VEB: ventricular ectopic, F: fusion, Q: unknown/paced)
BEAT_SYMBOLS = set("NLRBAaJSVrFejnE/fQ?")
AAMI = {
    **{s: "N" for s in "NLRBej"},
    **{s: "SVEB" for s in "AaJS"},
    **{s: "VEB" for s in "VE"},
    "F": "F",
    **{s: "Q" for s in "/fQ?"},
    "r": "VEB",   # R-on-T premature ventricular contraction
    "n": "SVEB",  # supraventricular escape beat
}

# The recordings have no absolute time, so all timestamps start from this fixed value
T0 = datetime(2026, 1, 1, 0, 0, 0, tzinfo=timezone.utc)


# ------------------------------------------------------------------ data loading
def records():
    """Names of the records listed in the database RECORDS file."""
    with open(os.path.join(DB, "RECORDS")) as f:
        return [r.strip() for r in f if r.strip()]


def load(rec):
    """Load one record.

    Returns the digital samples of lead MLII (or of the first channel when MLII is
    missing, i.e. records 102 and 104), the lead name, the converter gain and baseline,
    and the sample positions and symbols of the reference beat annotations."""
    r = wfdb.rdrecord(os.path.join(DB, rec), physical=False)
    names = r.sig_name
    ch = names.index("MLII") if "MLII" in names else 0
    sig = r.d_signal[:, ch].astype(np.int32)
    gain = float(r.adc_gain[ch])
    base = int(r.baseline[ch])
    a = wfdb.rdann(os.path.join(DB, rec), "atr")
    sym = np.array(a.symbol)
    smp = np.array(a.sample)
    is_beat = np.isin(sym, list(BEAT_SYMBOLS))   # drop rhythm and signal-quality annotations
    return sig, names[ch], gain, base, smp[is_beat], sym[is_beat]


def iso(t_seconds):
    """ISO 8601 timestamp (UTC) of a time offset in seconds from T0."""
    return (T0 + timedelta(seconds=float(t_seconds))).isoformat().replace("+00:00", "Z")


def dumps(obj):
    """Compact JSON, without indentation or spaces between elements."""
    return json.dumps(obj, separators=(",", ":"), ensure_ascii=False)


# ------------------------------------------------------------------ FHIR resource builders
MDC = "urn:oid:2.16.840.1.113883.6.24"                     # IEEE 11073 MDC code system
LEAD_CODE = {"MLII": ("131330", "MDC_ECG_ELEC_POTL_II")}   # lead II code, also used for MLII
CATEGORY = "http://terminology.hl7.org/CodeSystem/observation-category"
UCUM = "http://unitsofmeasure.org"


def ecg_observation(rec, lead, seg, start_s, gain, base):
    """Observation carrying an ECG segment in valueSampledData.

    The receiver computes value (mV) = origin + factor * data. With origin = -base/gain
    and factor = 1/gain, the integer samples of the converter are sent unchanged, so the
    mapping is lossless (eq. (1) of the article)."""
    code = [{"system": MDC, "code": "131328", "display": "MDC_ECG_ELEC_POTL"}]
    if lead in LEAD_CODE:
        c, d = LEAD_CODE[lead]
        code.append({"system": MDC, "code": c, "display": d})
    return {
        "resourceType": "Observation",
        "status": "final",
        "category": [{"coding": [{"system": CATEGORY, "code": "procedure"}]}],
        "code": {"coding": code, "text": f"ECG lead {lead}"},
        "subject": {"reference": f"Patient/mitdb-{rec}"},
        "device": {"reference": f"Device/gw-{rec}"},
        "effectiveDateTime": iso(start_s),
        "valueSampledData": {
            "origin": {"value": round(-base / gain, 6), "unit": "mV", "system": UCUM, "code": "mV"},
            "period": round(1000.0 / FS, 6),   # time between samples (ms)
            "factor": round(1.0 / gain, 6),
            "dimensions": 1,
            "data": " ".join(map(str, seg.tolist())),
        },
    }


def hr_observation(rec, start_s, dur_s, bpm):
    """Heart rate Observation (LOINC 8867-4, vital-signs profile) for one window."""
    return {
        "resourceType": "Observation",
        "meta": {"profile": ["http://hl7.org/fhir/StructureDefinition/vitalsigns"]},
        "status": "final",
        "category": [{"coding": [{"system": CATEGORY, "code": "vital-signs", "display": "Vital Signs"}]}],
        "code": {"coding": [{"system": "http://loinc.org", "code": "8867-4", "display": "Heart rate"}],
                 "text": "Heart rate"},
        "subject": {"reference": f"Patient/mitdb-{rec}"},
        "device": {"reference": f"Device/gw-{rec}"},
        "effectivePeriod": {"start": iso(start_s), "end": iso(start_s + dur_s)},
        "valueQuantity": {"value": round(float(bpm), 1), "unit": "beats/minute", "system": UCUM, "code": "/min"},
    }


def flag_interpretation(obs):
    """Add the Abnormal interpretation code to an ECG segment sent by the hybrid mode,
    so that the receiving system can forward it for review."""
    obs["interpretation"] = [{"coding": [{
        "system": "http://terminology.hl7.org/CodeSystem/v3-ObservationInterpretation",
        "code": "A", "display": "Abnormal"}],
        "text": "Irregular RR interval detected at the edge"}]
    return obs


def bundle(entries, key):
    """Transaction Bundle that creates the given resources.

    fullUrl is a UUID derived from `key`, so the messages, and therefore their sizes,
    are identical in every run."""
    return {
        "resourceType": "Bundle",
        "type": "transaction",
        "entry": [{"fullUrl": "urn:uuid:" + str(uuid.uuid5(uuid.NAMESPACE_URL, f"{key}/{i}")),
                   "resource": e,
                   "request": {"method": "POST", "url": "Observation"}}
                  for i, e in enumerate(entries)],
    }
