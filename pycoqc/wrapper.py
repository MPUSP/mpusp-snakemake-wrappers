__author__ = "Michael Jahn"
__copyright__ = "Copyright 2025, Michael Jahn"
__email__ = "jahn@mpusp.mpg.de"
__license__ = "MIT"


import os
from snakemake.shell import shell


# NOTE: Snakemake >=9 may inject modules that are not import-compatible with
# very old job Python runtimes (like pycoQC's historical py37 env). Accessing
# snakemake.params can trigger that import path, so we keep a compatibility
# fallback that reads the raw serialized params store.
# Fix was suggested by Github Copilot and implemented 2026-04-17.
def _get_param(name, default=""):
    try:
        return snakemake.params.get(name, default)
    except TypeError as e:
        if "unsupported operand type(s) for |" not in str(e):
            raise

        store = getattr(snakemake, "_params_store", None)
        if store is None:
            return default

        if isinstance(store, dict):
            return store.get(name, default)

        # Namedlist-like stores can expose names as attributes.
        if hasattr(store, name):
            return getattr(store, name)

        # Fallback for tuple-style key/value stores.
        for item in store:
            if isinstance(item, tuple) and len(item) == 2 and item[0] == name:
                return item[1]

        return default


# define default options for this wrapper
extra = _get_param("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# get input files
summary = snakemake.input.get("summary")
bam = snakemake.input.get("bam", "")
barcode = snakemake.input.get("barcode", "")


# prepare input
def concatenate_input(prefix, files):
    if files:
        files = " ".join(files) if isinstance(files, list) else files
        return f"--{prefix} " + files
    else:
        if prefix == "summary_file":
            raise ValueError("'summary' input is required for pycoQC.")
        return ""


summary = concatenate_input("summary_file", summary)
bam = concatenate_input("bam_file", bam)
barcode = concatenate_input("barcode_file", barcode)

# prepare output
report_html = snakemake.output.get("report_html")
report_json = snakemake.output.get("report_json")

# run pycioQC
shell(
    "pycoQC"
    " {summary}"  # required summary input
    " {bam}"  # optional bam input
    " {barcode}"  # optional barcode input
    " {extra}"  # additional options
    " -o {report_html}"  # html output
    " -j {report_json}"  # json output
    " {log}"
)
