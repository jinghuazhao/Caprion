import os

RWD = os.getcwd()

rule all:
    input:
        "results/tests.html"

rule report:
    input:
        markdown="scripts/tests.Rmd"
    output:
        outfile="results/tests.html"
    params:
        rwd = RWD,
        output_dir="results",
        output_name="tests.html"
    script:
        "../scripts/report.R"
