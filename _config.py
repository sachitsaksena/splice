import sys, os

PRJ_DIR = '/cluster/bh0085/prj/skipping4/'
SRC_DIR = PRJ_DIR + 'src/'

#######################################################
# Note: Directories should end in / always
#######################################################
DATA_DIR =os.path.join( PRJ_DIR + 'data/')
READS_DIR ="/cluster/bh0085/shortreads/190test"
OUT_PLACE =os.path.join( PRJ_DIR + 'out/')
RESULTS_PLACE = os.path.join(PRJ_DIR + 'results/')
QSUBS_DIR = os.path.join(PRJ_DIR + 'qsubs/')
#######################################################
#######################################################


EXP_DESIGN = os.path.join(DATA_DIR, "061318_exonskipping_library.csv")

#CONSTANTS
#SEQUENCING_DEFS = [e.split("\t") for e in """GEN00140700	120618 p2L + p2T SplAccLib preCas9 gDNA
#GEN00140701	121418 p2L + p2T SplAccLib +Cas9 rep1 gDNA
#GEN00140702	121418 p2L + p2T SplAccLib +Cas9 rep2 gDNA
#GEN00140703	120618 p2L + p2T SplAccLib preCas9 spliced RNA
#GEN00140704	121418 p2L + p2T SplAccLib +Cas9 rep1 spliced RNA
#GEN00140705	121418 p2L + p2T SplAccLib +Cas9 rep2 spliced RNA
#GEN00140706	120618 p2L + p2T SplAccLib preCas9 gDNA dictionary""".splitlines()]

FILE_FOLDERS_MESC = {
    "precas_unspliced" : "554",
    "unspliced1":'555',
    "unspliced2":'556',
    "precas_spliced":'557',
    "spliced1":'558',
    "spliced2":'559',
    "precas_gdna":'560',
    "gdna1":'561',
    "gdna2":'562',
}


FILE_FOLDERS_U2OS = {
    "precas_unspliced" : "563",
    "unspliced1":'564',
    "unspliced2":'565',
    "precas_spliced":'566',
    "spliced1":'567',
    "spliced2":'568',
    "precas_gdna":'569',
    "gdna1":'570',
    "gdna2":'571',
}
