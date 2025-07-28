__author__ = "Michael Jahn"
__copyright__ = "Copyright 2025, Michael Jahn"
__email__ = "jahn@mpusp.mpg.de"
__license__ = "MIT"


from os import path
from BCBio.GFF import GFFExaminer
from BCBio import GFF
from subprocess import getoutput
from snakemake.shell import shell

# define default options for this wrapper
threads = snakemake.threads
extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)
output_fasta = snakemake.output.get("fasta")
output_gff = snakemake.output.get("gff")
output_fai = snakemake.output.get("fai")
output_dir = path.commonpath([output_fasta, output_gff, output_fai])
if not output_dir:
    raise IOError("All output files should have the same non-empty base dir")

# get input files/params
default_sources = [
    {"RefSeq": "gene"},
    {"RefSeq": "pseudogene"},
    {"RefSeq": "CDS"},
    {"Protein Homology": "CDS"},
]
db = snakemake.params.get("database", "ncbi")
assembly = snakemake.params.get("assembly", "")
fasta = snakemake.input.get("fasta", "")
gff = snakemake.input.get("gff", "")
gff_source_types = snakemake.params.get("gff_source_types", default_sources)
messages = []


def check_fasta(fasta, messages=[]):
    with open(fasta, "r") as fasta_file:
        fasta = fasta_file.read()
    n_items = fasta.count(">")
    if n_items:
        messages += [f"Supplied fasta file '{fasta}' was found"]
        messages += [f"Supplied fasta file contains {n_items} items"]
    else:
        raise ValueError(
            "The supplied fasta file contains no valid entries starting with '>'"
        )
    return fasta, messages


def check_gff(input_gff, messages=[]):
    with open(input_gff, "r") as gff_file:
        gff_examiner = GFFExaminer()
        messages += [f"Supplied GFF file '{input_gff}' was found"]
        gff_summary = gff_examiner.available_limits(gff_file)
        messages += [
            f"Supplied GFF file contains the following items:",
            "--------------------",
        ]
        for item in gff_summary["gff_source_type"]:
            messages += [
                "-".join(item) + " : " + str(gff_summary["gff_source_type"][item])
            ]
    with open(input_gff, "r") as gff_file:
        new_gff = []
        gff_source_type = []
        for i in gff_source_types:
            gff_source_type += list(i.items())
        limits = dict(gff_source_type=gff_source_type)
        for rec in GFF.parse(gff_file, limit_info=limits):
            for recfeat in rec.features:
                rec_keys = recfeat.qualifiers.keys()
                if not "Name" in rec_keys:
                    if "locus_tag" in rec_keys:
                        recfeat.qualifiers["Name"] = recfeat.qualifiers["locus_tag"]
                    else:
                        raise ValueError(
                            "required fields 'Name','locus_tag' missing in *.gff file"
                        )
                else:
                    if "locus_tag" in rec_keys:
                        recfeat.qualifiers["trivial_name"] = recfeat.qualifiers["Name"]
                        recfeat.qualifiers["Name"] = recfeat.qualifiers["locus_tag"]
                if not "ID" in rec_keys:
                    if "locus_tag" in rec_keys:
                        recfeat.qualifiers["ID"] = recfeat.qualifiers["locus_tag"]
                    elif "Name" in rec_keys:
                        recfeat.qualifiers["ID"] = recfeat.qualifiers["Name"]
                    else:
                        raise ValueError(
                            "required fields 'ID','locus_tag' missing in *.gff file"
                        )
            new_gff += [rec]
    return new_gff, messages


if db.lower() == "ncbi":
    # Step 1: Query NCBI database
    shell(
        "(datasets summary genome accession "
        "{assembly} --as-json-lines | "
        "dataformat tsv genome "
        "--fields accession,annotinfo-release-date,organism-name > "
        "{output_dir}/ncbi_datasets_query.txt) "
        "{log}"
    )
    with open(path.join(output_dir, "ncbi_datasets_query.txt"), "r") as query:
        ncbi_result = query.read()
    if ncbi_result.startswith("Error"):
        raise ValueError(
            "The supplied refseq/genbank ID was not valid. "
            + "Example for correct input: 'GCF_000009045.1'"
        )
    if len(ncbi_result) == 0:
        raise ValueError(
            "The result from fetching NCBI genome data has zero length. "
            + "Please check your internet connection!"
        )

    ncbi_genome = [
        i.split("\t")
        for i in ncbi_result.split("\n")
        if not (i.startswith("New version") or i.startswith("Warning"))
    ]
    ncbi_genome = dict(zip(ncbi_genome[0], ncbi_genome[1]))
    messages += ["Found the following genome(s):"]
    for k in ncbi_genome.keys():
        messages += ["{0}: {1}".format(k, ncbi_genome.get(k))]
    refseq_id = ncbi_genome.get("Assembly Accession")
    if not refseq_id.startswith("GCF_"):
        raise ValueError(f"The RefSeq ID '{refseq_id}' has no valid format")

    # Step 2: Retrieve dataset from NCBI
    shell(
        "(datasets download genome accession "
        "{refseq_id} "
        "--filename {output_dir}/database.zip "
        "--include genome,gff3; "
        "cd {output_dir}; "
        "unzip -o database.zip; "
        "rm database.zip) "
        "{log}"
    )

    # Step 3: Copy files to final output location
    shell(
        "cp {output_dir}/ncbi_dataset/data/{refseq_id}/*.fna {output_fasta} {log}; "
        "cp {output_dir}/ncbi_dataset/data/{refseq_id}/genomic.gff {output_gff} {log}"
    )

    # Step 4: Index FASTA file
    shell("samtools faidx {output_fasta} {log}")

    # check and export files
    fasta, messages = check_fasta(output_fasta, messages)
    gff, messages = check_gff(output_gff, messages)
    with open(output_gff, "w") as gff_out:
        GFF.write(gff, gff_out)

elif db.lower() == "manual":
    if not path.exists(fasta):
        raise IOError("The parameter 'fasta' is not a valid path to a FASTA file")
    elif not path.exists(gff):
        raise IOError("The parameter 'gff' is not a valid path to a GFF/GTF file")
    else:
        # import and check files
        fasta, messages = check_fasta(fasta, messages)
        gff, messages = check_gff(gff, messages)
        # export fasta and gff files
        with open(output_fasta, "w") as fasta_out:
            fasta_out.write(fasta)
        with open(output_gff, "w") as gff_out:
            GFF.write(gff, gff_out)
        # index FASTA file
        shell("samtools faidx {output_fasta} {log}")
else:
    raise ValueError("The parameter 'database' is none of 'ncbi', 'manual'")
