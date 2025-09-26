import yaml

# Load config
configfile: "config.yaml"
QUNCY_SRC_ID   = config["quincy_source"]["gdrive_id"]
QUNCY_SRC_DIR  = config["quincy_source"]["local_dir"]

# Data info from config
DATA_URL   = config["data"]["gdrive_url"]
DATA_DIR   = config["data"]["local_dir"]

# Final workflow target
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
        directory("results/python_test_imgs")
 


# Rule to build the Fortran project
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

# Rule to prepare the Python environment
rule python_env:
    output: touch("env/.installed")
    conda: "env/environment.yaml"

# Fortran test
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

# Reporting rule: print results and mark completion
rule report_tests:
    input:
        "tests/fortran_works.txt",
    output:
        touch("tests/.report_done")
    run:
        for f in input:
            print(f"\n--- {f} ---")
            with open(f) as fh:
                print(fh.read())


# Rule to fetch data from Google Drive
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

        # download from Google Drive
        gdown --id {QUNCY_SRC_ID} -O $tmpzip

        # clean old directory
        rm -rf {output}
        mkdir -p {output}

        # extract and strip the top-level folder
        unzip -q $tmpzip -d {output}_tmp
        mv {output}_tmp/*/* {output}/ || mv {output}_tmp/* {output}/
        rm -rf {output}_tmp
        rm -f $tmpzip
        """



rule record_python:
    input:
        "env/.installed"
    output:
        "results/python_path.txt"
    conda:
        "env/environment.yaml"
    shell:
        """
        echo $CONDA_PREFIX/bin/python > {output}
        """


rule record_data:
    input:
        rules.get_data.output
    output:
        "results/data_files.txt"
    shell:
        """
        # Only record the two expected files if they exist
        realpath {input}/forcing/transient/ATTO_t_1901-2023.dat > {output}
        realpath {input}/forcing/static/ATTO_s_2000-2023.dat >> {output}
        realpath {input}/obs.csv >> {output}
        """

rule record_quincy:
    input:
        "build/quincy_run"
    output:
        "results/quincy_path.txt"
    shell:
        """
        realpath {input} > {output}
        """


rule test_python_script:
    input:
        "tests/.report_done",        # after Fortran test + report
        "env/.installed",            # environment ready
        "build/quincy_run",          # Fortran build done
        directory(DATA_DIR),         # data available
        ".submodules/QPy_checked"    # submodule checked out
    output:
        directory("results/python_test_imgs")
    conda:
        "env/environment.yaml"
    shell:
        """
        mkdir -p results/python_test_imgs
        $CONDA_PREFIX/bin/python -u QPy/science/atto_excercises/00_run_static_quincy.py

        # Ensure at least one PDF exists, otherwise fail
        test $(ls results/python_test_imgs/*.pdf | wc -l) -gt 0
        """

rule update_submodules:
    output: ".submodules/QPy_checked"
    shell:
        r"""
        set -euo pipefail

        # Make sure local config matches .gitmodules
        git submodule sync --recursive

        # Ensure the submodule worktree exists locally
        git submodule update --init --recursive QPy

        # Ensure the desired branch is checked out & up to date
        git -C QPy fetch origin feature/atto_summer_phyd
        git -C QPy checkout -B feature/atto_summer_phyd origin/feature/atto_summer_phyd
        git -C QPy pull --ff-only || true

        mkdir -p .submodules
        touch {output}
        """

