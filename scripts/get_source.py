import sys
import zipfile
import tempfile
import shutil
from pathlib import Path
import gdown

# Args from Snakemake
gdrive_id = snakemake.params.gdrive_id
output_dir = Path(snakemake.output[0])
log_file = Path(snakemake.log[0])

def log(msg):
    with log_file.open("a") as f:
        f.write(f"[INFO] {msg}\n")
    print(f"[INFO] {msg}", flush=True)

def main():
    log(f"Downloading Google Drive file ID={gdrive_id}")

    # Temporary zip file
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpzip = Path(tmpdir) / "source.zip"

        # Download
        url = f"https://drive.google.com/uc?id={gdrive_id}"
        gdown.download(url, str(tmpzip), quiet=False)

        log(f"Downloaded to {tmpzip}, size={tmpzip.stat().st_size} bytes")

        # Check if valid zip
        if not zipfile.is_zipfile(tmpzip):
            raise RuntimeError(f"Downloaded file is not a valid zip: {tmpzip}")

        # Clean old output
        if output_dir.exists():
            shutil.rmtree(output_dir)
        output_dir.mkdir(parents=True)

        # Extract into temp dir
        tmp_extract = Path(tmpdir) / "extract"
        with zipfile.ZipFile(tmpzip, "r") as zf:
            zf.extractall(tmp_extract)

        log(f"Extracted contents: {[p.name for p in tmp_extract.iterdir()]}")

        # Move files into output (handles both flat and nested zips)
        # If there is exactly one top-level folder, flatten it
        items = list(tmp_extract.iterdir())
        if len(items) == 1 and items[0].is_dir():
            src = items[0]
        else:
            src = tmp_extract

        for item in src.iterdir():
            dest = output_dir / item.name
            if item.is_dir():
                shutil.copytree(item, dest)
            else:
                shutil.copy2(item, dest)

        log(f"Final contents of {output_dir}: {[p.name for p in output_dir.iterdir()]}")

if __name__ == "__main__":
    main()
