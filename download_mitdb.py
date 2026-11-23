"""Download the MIT-BIH Arrhythmia Database v1.0.0 from PhysioNet into data/.

Source: https://physionet.org/content/mitdb/1.0.0/ (https://doi.org/10.13026/C2F305),
Open Data Commons Attribution License v1.0. The integrity of the files is checked
against the archive's own SHA256SUMS.txt by 05_checks.py."""
import os
import urllib.request
import zipfile

URL = "https://physionet.org/static/published-projects/mitdb/mit-bih-arrhythmia-database-1.0.0.zip"
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
ZIP = os.path.join(DATA, "mit-bih-arrhythmia-database-1.0.0.zip")


def main():
    os.makedirs(DATA, exist_ok=True)
    if not os.path.exists(ZIP):
        print("downloading", URL)
        urllib.request.urlretrieve(URL, ZIP)
    with zipfile.ZipFile(ZIP) as z:
        z.extractall(DATA)
    print("records in", os.path.join(DATA, "mit-bih-arrhythmia-database-1.0.0"))


if __name__ == "__main__":
    main()
