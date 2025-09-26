import yaml

# Load config
configfile: "config.yaml"
CMAKE_SOURCE = config["cmake_source"]
# Repository info from config
REPO_URL   = config["repo"]["url"]
BRANCH     = config["repo"]["branch"]
CLONE_DIR  = config["repo"]["clone_dir"]

# Data info from config
DATA_URL   = config["data"]["gdrive_url"]
DATA_DIR   = config["data"]["local_dir"]

# Final workflow target
rule all:
    input:
        "build/quincy_run",
        "env/.installed",
        "tests/fortran_works.txt",
        "tests/.report_done", 
        directory(CLONE_DIR),
        directory(DATA_DIR),
        "results/python_path.txt",
        "results/quincy_path.txt",
        "results/data_files.txt",
        directory("results/python_test_imgs")
 


# Rule to build the Fortran project
rule build_fortran:
    output: "build/quincy_run"
    conda: "env/environment.yaml"
    shell:
        """
        cmake -S {CMAKE_SOURCE} -B build \
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


rule clone_repo:
    output:
        directory(CLONE_DIR)
    shell:
        """
        if [ -d "{output}" ]; then
            cd {output} && git fetch && git checkout {BRANCH} && git pull
        else
            git clone --branch {BRANCH} {REPO_URL} {output}
        fi
        """


# Rule to fetch data from Google Drive
rule get_data:
    output: directory(DATA_DIR)
    conda: "env/environment.yaml"
    shell:
        """
        mkdir -p {output}
        gdown --folder {DATA_URL} -O {output}
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
        "tests/.report_done",        # only run after Fortran test + report finished
        "env/.installed",            # environment must be ready
        "build/quincy_run",          # ensures Fortran build succeeded
        directory(DATA_DIR)          # ensures data is available
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


