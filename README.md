# Electrocardiogram encoding in HL7 FHIR for remote monitoring: experiment code

Code for the article "Electrocardiogram encoding in HL7 FHIR for remote monitoring: data volume,
server cost and edge accuracy" (M. Ileana, R. Sfat, C.V. Marian), submitted to
Rev. Roum. Sci. Techn. – Électrotechn. et Énerg. Running the scripts in order reproduces every number,
table and figure of the article.

## What the experiment does

A single-lead ECG from the MIT-BIH Arrhythmia Database is processed on an emulated edge gateway, encoded
into HL7 FHIR resources and sent over HTTP/1.1 to a FHIR endpoint that validates every message.
Six encodings are compared:

| ID  | Content and message unit |
|-----|--------------------------|
| S1  | raw ECG, one Observation per 1 s, one request per Observation |
| S1b | raw ECG, one Observation per 10 s, one request per Observation |
| S2  | raw ECG, transaction Bundle of six 10-s Observations per 60 s |
| S3  | as S2, gzip |
| S4  | heart rate only, Bundle of 10-s HR Observations per 60 s, gzip |
| S5  | S4 plus the 10-s ECG segments that contain a premature beat (θ = 0.85), gzip |

S1 → S1b changes only the length of the Observation; S1b → S2 changes only the grouping in Bundles.

## Data

MIT-BIH Arrhythmia Database v1.0.0, PhysioNet: https://physionet.org/content/mitdb/1.0.0/
(https://doi.org/10.13026/C2F305), Open Data Commons Attribution License v1.0.
The data are not included in this repository. Download them with

    python download_mitdb.py

which places the records in `data/mit-bih-arrhythmia-database-1.0.0/`. Another location can be given
with the `MITDB` environment variable.

## Environment

Python 3.13 and the packages in `requirements.txt`:

    pip install -r requirements.txt

Host used for the article: Intel Xeon 2.10 GHz, two virtual cores, 8 GB of memory.
Figures use the TeX Gyre Termes font if it is installed (folder given by `TEXGYRE_DIR` or the usual TeX
paths), otherwise the default serif font.

`05_checks.py` compares record 100 with the copy distributed with wfdb-python. Clone
https://github.com/MIT-LCP/wfdb-python into `wfdbrepo/` or set `WFDB_SAMPLE` to its `sample-data` folder.

## Running

    python download_mitdb.py
    python 01_detect.py
    python 02_encode_send.py
    python 03_extra.py
    python 04_figures.py
    python 05_checks.py

| Script | What it does | Output |
|---|---|---|
| `common.py` | data loading, FHIR resource builders, AAMI beat classes | – |
| `01_detect.py` | streaming Pan-Tompkins QRS detection (60-s blocks), scored against the annotations (150 ms) | `results/detection.json`, `results/peaks_*.npy` |
| `02_encode_send.py` | builds S1–S5, sends them to a local validating FHIR endpoint; data volume, timing, HR error, ectopic beat capture, threshold sweep, lossless check | `results/encode_send.json` |
| `03_extra.py` | per-record S5 volume, share of windows with ectopic beats, false alarms without beat annotations, whole-record detection | `results/extra.json` |
| `04_figures.py` | Fig. 1, 3 and 4 and the Fig. 4 regression statistics | `fig/*.png`, `results/fig4_stats.json` |
| `05_checks.py` | validation negative control, lead per record, SHA-256 check of the database, gain and baseline of every record, SampledData parameters, record 100 against the wfdb copy | `results/checks.json` |

The `results/` folder holds the outputs of the run reported in the article (strategies S1–S5);
running the scripts again overwrites them and adds S1b.

Data volumes, request counts, accuracy values and figures are deterministic for a given zlib
library. The size of the gzip-compressed messages (S3–S5) depends on the zlib implementation: with
zlib-ng, used by some recent Python builds, the compressed messages are a few percent smaller than with
the reference zlib. The version used is written to `results/encode_send.json` (`zlib`). Timing values (round-trip
time, server time, gateway processor time) depend on the machine and vary between runs.
`results/rerun_comparison.json` records repeated runs on the same type of virtual machine: the
deterministic outputs were identical, and the timing values differed by up to 31% for per-hour totals,
36% for medians and 46% for 95th percentiles.

## Notes

- The gateway and the endpoint run as separate processes on one host and communicate only over HTTP
  (loopback); a single client sends the requests one after another. Times do not include network delay
  or TLS.
- Validation is structural (fhir.resources 8.3.0, FHIR R4B models); terminology bindings are not checked.
- Timestamps are synthetic (origin 2026-01-01T00:00:00Z), because the recordings have no absolute time.

## License

Code: MIT (see `LICENSE`). Data: see the PhysioNet license above.
