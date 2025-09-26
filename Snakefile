import yaml
import json
from pathlib import Path

# -------------------------------------------------------
# Helper: update VS Code interpreter from a path file
# -------------------------------------------------------
def update_vscode(python_path_file):
    from pathlib import Path
    import json

    with open(python_path_file) as f:
        python_path = f.read().strip()

    # explicit submodule path
    settings_file = Path("QPy/.vscode/settings.json")

    # load existing settings if present
    settings = {}
    if settings_file.exists():
        try:
            with open(settings_file) as f:
                settings = json.load(f)
        except json.JSONDecodeError:
            print(f"[Snakemake] Warning: {settings_file} was invalid, starting fresh.")

    # update/add interpreter key
    settings["python.defaultInterpreterPath"] = python_path

    # ensure directory exists
    settings_file.parent.mkdir(exist_ok=True)
    with open(settings_file, "w") as f:
        json.dump(settings, f, indent=4)

    print(f"[Snakemake] Wrote interpreter to {settings_file}")
    print(f"[Snakemake] VS Code interpreter set to {python_path}")




# -------------------------------------------------------
# Load config
# -------------------------------------------------------
configfile: "config.yaml"
QUNCY_SRC_ID   = config["quincy_source"]["gdrive_id"]
QUNCY_SRC_DIR  = config["quincy_source"]["local_dir"]

DATA_URL   = config["data"]["gdrive_url"]
DATA_DIR   = config["data"]["local_dir"]

# -------------------------------------------------------
# Final workflow target
# -------------------------------------------------------
rule all:
    input:
        directory(QUNCY_SRC_DIR),
        "build/quincy_run",
        ".submodules/QPy_checked",
        "env/.installed",
        "tests/fortran_works.txt",
        "tests/.report_done",
        directory(DATA_DIR),
        "results/python_path.txt",
        "results/quincy_path.txt",
        "results/data_files.txt",
        directory("results/python_test_imgs"),
        "results/.vscode_synced"   # <--- new sync marker

# -------------------------------------------------------
# Rules
# -------------------------------------------------------

rule build_fortran:
    input: directory(QUNCY_SRC_DIR)
    output: "build/quincy_run"
    conda: "env/environment.yaml"
    shell:
        """
        cmake -S {input} -B build \
              -DCMAKE_Fortran_COMPILER=gfortran \
              -DCMAKE_PREFIX_PATH=$CONDA_PREFIX
        cmake --build build --parallel
        """

rule python_env:
    output: touch("env/.installed")
    conda: "env/environment.yaml"

rule test_fortran:
    output: "tests/fortran_works.txt"
    conda: "env/environment.yaml"
    shell:
        """
        mkdir -p tests
        echo 'program test; print *, "Fortran works!"; end program test' > tests/test.f90
        gfortran tests/test.f90 -o tests/test_prog
        ./tests/test_prog > {output}
        echo "✅ Fortran test completed" >> {output}
        """

rule report_tests:
    input: "tests/fortran_works.txt"
    output: touch("tests/.report_done")
    run:
        for f in input:
            print(f"\n--- {f} ---")
            with open(f) as fh:
                print(fh.read())

rule get_data:
    output: directory(DATA_DIR)
    conda: "env/environment.yaml"
    shell:
        """
        mkdir -p {output}
        gdown --folder {DATA_URL} -O {output}
        """

rule get_source:
    output: directory(QUNCY_SRC_DIR)
    conda: "env/environment.yaml"
    shell:
        """
        tmpzip=$(mktemp /tmp/sourceXXXX.zip)
        gdown --id {QUNCY_SRC_ID} -O $tmpzip
        rm -rf {output}
        mkdir -p {output}
        unzip -q $tmpzip -d {output}_tmp
        mv {output}_tmp/*/* {output}/ || mv {output}_tmp/* {output}/
        rm -rf {output}_tmp
        rm -f $tmpzip
        """

rule record_python:
    input: "env/.installed"
    output: "results/python_path.txt"
    conda: "env/environment.yaml"
    shell:
        """
        echo $CONDA_PREFIX/bin/python > {output}
        """

rule record_data:
    input: rules.get_data.output
    output: "results/data_files.txt"
    shell:
        """
        realpath {input}/forcing/transient/ATTO_t_1901-2023.dat > {output}
        realpath {input}/forcing/static/ATTO_s_2000-2023.dat >> {output}
        realpath {input}/obs.csv >> {output}
        """

rule record_quincy:
    input: "build/quincy_run"
    output: "results/quincy_path.txt"
    shell:
        "realpath {input} > {output}"

rule test_python_script:
    input:
        "tests/.report_done",
        "env/.installed",
        "build/quincy_run",
        directory(DATA_DIR),
        ".submodules/QPy_checked"
    output: directory("results/python_test_imgs")
    conda: "env/environment.yaml"
    shell:
        """
        mkdir -p results/python_test_imgs
        $CONDA_PREFIX/bin/python -u QPy/science/atto_excercises/00_run_static_quincy.py
        test $(ls results/python_test_imgs/*.pdf | wc -l) -gt 0
        """

rule update_submodules:
    output: ".submodules/QPy_checked"
    shell:
        r"""
        set -euo pipefail
        git submodule sync --recursive
        git submodule update --init --recursive QPy
        git -C QPy fetch origin feature/atto_summer_phyd
        git -C QPy checkout -B feature/atto_summer_phyd origin/feature/atto_summer_phyd
        git -C QPy pull --ff-only || true
        mkdir -p .submodules
        touch {output}
        """

# -------------------------------------------------------
# New rule: sync VS Code
# -------------------------------------------------------
rule sync_vscode:
    input: "results/python_path.txt"
    output: touch("results/.vscode_synced")
    run:
        update_vscode(input[0])
