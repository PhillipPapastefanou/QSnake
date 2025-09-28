# QSnake

QSnake is a Snakemake-based workflow for building and running the **Quincy** model together with Python tooling from the **QPy** submodule.  

---

## 📋 Prerequisites

Before running the workflow, ensure you have:

- **[Conda / Mamba](https://docs.conda.io/en/latest/miniconda.html)** (recommended for environment management)
- **[Snakemake](https://snakemake.readthedocs.io/)**  
  Install via Conda:
  ```bash
  conda create -n snakemake -c conda-forge -c bioconda snakemake
  conda activate snakemake
  ```
- **CMake** (for building Fortran code)
- **Fortran compiler** (`gfortran` recommended)
- **Git** (for cloning and submodules)

💡 On **Windows**, the easiest setup is via **WSL2 (Ubuntu)** so you have access to `gfortran` and other Unix tools.

---

## 🚀 General Setup

1. **Clone the repository**
   ```bash
   git clone https://github.com/PhillipPapastefanou/QSnake
   cd QSnake
   ```

2. **Initialize and update the QPy submodule**
   ```bash
   git submodule update --init --recursive
   ```

3. **Run the build with Snakemake**
   ```bash
   snakemake --use-conda --cores 2 --latency-wait 10
   ```

---

## 🛠 How it works

- **Snakemake** orchestrates the workflow, creating environments as needed.  
- **QPy** is included as a submodule and kept up-to-date by the workflow.  
- **Fortran sources** (Quincy) are built with CMake and gfortran.  
- **Python scripts** and tests ensure the build and environment are working correctly.  

---

## 📂 Project Structure

```
QSnake/
├── Snakefile          # Main workflow
├── config.yaml        # Configuration (Google Drive IDs, paths, etc.)
├── env/               # Conda environment files
├── scripts/           # Helper Python scripts
├── results/           # Outputs (build paths, tests, plots)
└── QPy/               # Python submodule
```

---

## ✅ Testing

After the build completes, you can verify that everything works:

```bash
snakemake tests/.report_done --use-conda
```

This runs a minimal Fortran compilation and Python checks.

---

## 📖 Notes

- Make sure your Google Drive links in `config.yaml` are accessible (set to “Anyone with the link → Viewer”).
- On Windows, WSL2 is strongly recommended for compatibility.
- Logs from Snakemake rules are stored under `logs/`.

---

## 📜 License

See [LICENSE.md](LICENSE.md) for details.
