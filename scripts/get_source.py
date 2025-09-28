import zipfile
import tempfile
import shutil
from pathlib import Path
import gdown

# Snakemake injects these
gdrive_id = snakemake.params.gdrive_id
output_dir = Path(snakemake.output[0])
log_file = Path(snakemake.log[0])

def log(msg: str):
    """Write to log file and stdout."""
    with log_file.open("a") as f:
        f.write(f"[INFO] {msg}\n")
    print(f"[INFO] {msg}", flush=True)

def safe_extract(zf: zipfile.ZipFile, target: Path):
    """
    Extracts a zipfile while restoring permissions (executable bits, symlinks).
    """
    for member in zf.infolist():
        target_path = target / member.filename

        if member.is_dir():
            target_path.mkdir(parents=True, exist_ok=True)
            continue

        # Extract the file
        zf.extract(member, target)

        # Restore original permissions if present
        mode = member.external_attr >> 16
        if mode:
            target_path.chmod(mode)

        # Handle symlinks (zip stores symlinks as special files with mode 0o120000)
        if (mode & 0o170000) == 0o120000:
            link_target = target_path.read_text()
            target_path.unlink()
            target_path.parent.mkdir(parents=True, exist_ok=True)
            target_path.symlink_to(link_target)

def main():
    log(f"Downloading Google Drive file ID={gdrive_id}")

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpzip = Path(tmpdir) / "source.zip"

        # Download file
        url = f"https://drive.google.com/uc?id={gdrive_id}"
        gdown.download(url, str(tmpzip), quiet=False)

        log(f"Downloaded to {tmpzip}, size={tmpzip.stat().st_size} bytes")

        if not zipfile.is_zipfile(tmpzip):
            raise RuntimeError(f"Downloaded file is not a valid zip: {tmpzip}")

        # Clean old output
        if output_dir.exists():
            shutil.rmtree(output_dir)
        output_dir.mkdir(parents=True)

        # Extract to temp dir
        tmp_extract = Path(tmpdir) / "extract"
        tmp_extract.mkdir()
        with zipfile.ZipFile(tmpzip, "r") as zf:
            safe_extract(zf, tmp_extract)

        log(f"Top-level extracted entries: {[p.name for p in tmp_extract.iterdir()]}")

        # Flatten if the archive has a single top-level folder
        items = list(tmp_extract.iterdir())
        if len(items) == 1 and items[0].is_dir():
            src = items[0]
        else:
            src = tmp_extract

        # Copy into output_dir
        for item in src.iterdir():
            dest = output_dir / item.name
            if item.is_dir():
                shutil.copytree(item, dest)
            else:
                shutil.copy2(item, dest)

        # Sanity check for CMakeLists.txt
        cmake_file = output_dir / "CMakeLists.txt"
        if not cmake_file.exists():
            raise RuntimeError(
                f"Expected CMakeLists.txt at {cmake_file}, "
                f"but found {[p.name for p in output_dir.iterdir()]}"
            )

        log(f"Final contents of {output_dir}: {[p.name for p in output_dir.iterdir()]}")

if __name__ == "__main__":
    main()
