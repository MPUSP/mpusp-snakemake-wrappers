[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg)](https://snakemake.readthedocs.io)

# MPUSP Snakemake wrappers

## Purpose

This repository is a collection of custom [Snakemake](https://snakemake.readthedocs.io) wrappers. The Snakemake community maintains a wide variety of [Wrappers](https://github.com/snakemake/snakemake-wrappers) for recurring tasks in the life sciences.

A wrapper is simply a standardized script with a [pinned software environment](https://snakemake-wrappers.readthedocs.io/en/stable/contributing.html#conda-environment-for-development) that can be re-used in existing workflows. Wrappers often simply provide a reproducible interface to existing tools. How wrappers work is covered in detail on the [Snakemake Docs](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#wrappers) and the [Snakemake Wrapper Homepage](https://snakemake-wrappers.readthedocs.io/en/stable/contributing.html)

## Use wrappers in a project

Wrappers can be included in any workflow using the `wrapper` directive instead of `script` or `shell`.
This is an example for using the `GFFCOMPARE` wrapper. Note that a GitHub URL needs to be provided that points to the `raw.githubusercontent.com` destination of a wrapper, not the rendered HTML.

```bash
rule test_gffcompare:
    input:
        gff=["file1.gff", "file2.gff", "file3.gff"],
        # optional input files
        # reference="", # reference transcriptome in GFF/GTF format
        # fasta="", # reference genome sequence in FASTA format
    output:
        # it is recommended to use the `multiext` statement, as all output files are
        # created using a common prefix
        multiext("gffcmp", ".combined.gtf", ".loci", ".stats", ".tracking"),
    threads: 1
    log:
        "logs/gffcompare.log",
    params:
        # To combine features from all GFF input files without redundancy, use `-C`
        extra="-C",
    wrapper:
        "https://raw.githubusercontent.com/MPUSP/mpusp-snakemake-wrappers/refs/heads/main/gffcompare"
```

## Adding wrappers to this repository

Wrappers that are to be added to this repository shoul alreadz follow the standards of the official Snakemake wrapper repository. This way it will be easy to transfer them to the official wrapper collection if desired.

## License

- Wrappers in this repository are MIT-licensed
- Note: all wrapped/external software is licensed under its own terms
- Wrappers only provide an interface to call external software using the Snakemake workflow manager
