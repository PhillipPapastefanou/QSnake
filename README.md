# General setup

First run:
git clone https://github.com/PhillipPapastefanou/QSnake

then make sure that the QPy submodule is being pulled and up to date:
git submodule update --init --recursive

Finally, run the build script:  
snakemake --use-conda --cores 2 --latency-wait 10
