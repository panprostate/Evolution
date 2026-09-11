#!/usr/bin/env python
from cyvcf2 import VCF, Writer
from pathlib import Path
from multiprocessing import Pool

# Variables
cores = 10
snvdir = Path("data/processed/mutation_timer/snvs")
indeldir = Path("data/processed/mutation_timer/indels")

if not snvdir.exists() or not indeldir.exists():
    raise ValueError("SNV or INDEL dir does not exist")

test_vcf = Path("data/processed/mutation_timer/snvs/PPCG0002a_DNA_vs_PPCG0002b_DNA_mutation_timer.vcf")

# Set paths for clonal and subclonal dirs
sbs_dbs_prefix = "outputs/00_preprocessing/sigmatgen"
sbs_clonal_out = Path(sbs_dbs_prefix, "sbs/clonal")
sbs_subclonal_out = Path(sbs_dbs_prefix, "sbs/subclonal")
sbs_early_out = Path(sbs_dbs_prefix, "sbs/early")
sbs_late_out = Path(sbs_dbs_prefix, "sbs/late")
dbs_clonal_out = Path(sbs_dbs_prefix, "dbs/clonal")
dbs_subclonal_out = Path(sbs_dbs_prefix, "dbs/subclonal")
dbs_early_out = Path(sbs_dbs_prefix, "dbs/early")
dbs_late_out = Path(sbs_dbs_prefix, "dbs/late")
indel_clonal_out = Path(sbs_dbs_prefix, "id/clonal")
indel_subclonal_out = Path(sbs_dbs_prefix, "id/subclonal")
indel_early_out = Path(sbs_dbs_prefix, "id/early")
indel_late_out = Path(sbs_dbs_prefix, "id/late")
output_directories = [ sbs_clonal_out, sbs_subclonal_out, sbs_early_out, sbs_late_out, dbs_clonal_out,
                   dbs_subclonal_out, dbs_early_out, dbs_late_out, indel_clonal_out, indel_subclonal_out,
                    indel_early_out, indel_late_out ]

# Create outdirs if they don't exist
for output_dir in output_directories:
    if not output_dir.exists():
        output_dir.mkdir(parents=True, exist_ok=True)


# Setup functions ----
def is_dbs(var, next_var):
    if not (next_var is None):
        if var.CHROM == next_var.CHROM and var.start == (next_var.start - 1) and var.INFO.get("CLS") == next_var.INFO.get("CLS"):
            return True
    return False

def write_cls(variant, cfl, sfl, efl, lfl):
    # Write variant to clonal, subclonal, early or late files.
    # File handles must be Writer object instances
    cls = variant.INFO.get("CLS")
    if cls in ['clonal [NA]', 'clonal [early]', 'clonal [late]']:
        cfl.write_record(variant)
    if cls in ['subclonal']:
        sfl.write_record(variant)
    if cls in ['clonal [early]']:
        efl.write_record(variant)
    if cls in ['clonal [late]']:
        lfl.write_record(variant)
    

def parse_sbs_dbs_vcf(vcf_file):

    print("START_SBS: {}".format(vcf_file))

    # Loop variants and write to outfiles
    vcf = VCF(vcf_file)

    # Setup output files for the sample
    sbs_clf = Writer(Path(sbs_clonal_out, vcf_file.name), vcf)
    sbs_sfl = Writer(Path(sbs_subclonal_out, vcf_file.name), vcf)
    sbs_efl = Writer(Path(sbs_early_out, vcf_file.name), vcf)
    sbs_lfl = Writer(Path(sbs_late_out, vcf_file.name), vcf)
    dbs_clf = Writer(Path(dbs_clonal_out, vcf_file.name), vcf)
    dbs_sfl = Writer(Path(dbs_subclonal_out, vcf_file.name), vcf)
    dbs_efl = Writer(Path(dbs_early_out, vcf_file.name), vcf)
    dbs_lfl = Writer(Path(dbs_late_out, vcf_file.name), vcf)

    for variant in vcf:
        # WARNING: Expects all variants to be single bases

        # Get the next variant. Note that this consumes it.
        try:
            next_var = next(vcf)
        except StopIteration:
            next_var = None

        # If DBS, write to DBS file and continue
        dbs_call = is_dbs(variant, next_var)
        if(dbs_call):
            # Write to DBS file for cls
            write_cls(variant, dbs_clf, dbs_sfl, dbs_efl, dbs_lfl)
            write_cls(next_var, dbs_clf, dbs_sfl, dbs_efl, dbs_lfl)
            continue

        # Otherwise, write both to SBS files for class
        write_cls(variant, sbs_clf, sbs_sfl, sbs_efl, sbs_lfl)
        if not (next_var is None):
            write_cls(next_var, sbs_clf, sbs_sfl, sbs_efl, sbs_lfl)
    
    # Close all file handles
    sbs_clf.close(); sbs_sfl.close(); sbs_efl.close(); sbs_lfl.close()
    dbs_clf.close(); dbs_sfl.close(); dbs_efl.close(); dbs_lfl.close()
    vcf.close()

    print("END_SBS: {}".format(vcf_file))
    return(vcf_file.name)

def parse_indel_vcf(vcf_file):

    print("START_INDEL: {}".format(vcf_file))

    # Loop variants and write to outfiles
    vcf = VCF(vcf_file)

    # Setup output files for the sample
    indel_clf = Writer(Path(indel_clonal_out, vcf_file.name), vcf)
    indel_sfl = Writer(Path(indel_subclonal_out, vcf_file.name), vcf)
    indel_efl = Writer(Path(indel_early_out, vcf_file.name), vcf)
    indel_lfl = Writer(Path(indel_late_out, vcf_file.name), vcf)

    for variant in vcf:
        # Write to indel file for cls
        write_cls(variant, indel_clf, indel_sfl, indel_efl, indel_lfl)
    
    # Close all file handles
    indel_clf.close(); indel_sfl.close(); indel_efl.close(); indel_lfl.close()
    vcf.close()
    
    print("END_INDEL: {}".format(vcf_file))
    return(vcf_file.name)


# Run parser on all mutation_timer VCF files in parallel
snv_vcfs = list(Path(snvdir).glob("*.vcf"))
indel_vcfs = list(Path(indeldir).glob("*.vcf"))

with Pool(cores) as p:
    p.map(parse_sbs_dbs_vcf, snv_vcfs)

with Pool(cores) as p:
    p.map(parse_indel_vcf, indel_vcfs)


# Run SigProfilerMatrixGenerator on all output directories
# Setup:
# Install reference for SigProfilerMatrixGenerator
#from SigProfilerMatrixGenerator import install as genInstall
#genInstall.install('GRCh37', rsync=False, bash=True)

# Generate matrices
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

# Call on each output directory
for dir in output_directories:
    project = dir.parent.name + "::" + dir.name
    print("START_MATRIX: {}".format(project))
    dirpath = str(dir) + "/"
    matrices = matGen.SigProfilerMatrixGeneratorFunc(
        project, "GRCh37", dirpath,plot=False,
        exome=False, bed_file=None,
        chrom_based=False, tsb_stat=False,
        seqInfo=False, cushion=100)