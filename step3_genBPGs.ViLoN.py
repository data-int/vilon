#!/usr/bin/python3
# -*- coding: utf-8 -*-
# run with: snakemake -s step3_genBPGs.ViLoN.py -j 20
# create a DAG: snakemake -s step3_genBPGs.ViLoN.py --dag | dot -Tpdf > dag.pdf # visualization of this particular computation; does not always work (needs more dependencies; see Snakemake online, for more details)

from os import listdir
from os.path import isfile, join

## PROGZ ##
PNL='tools/genBPGs.pl'
PNL_MODE='abs'

## VARIABLES ##
CAN='vilon.online'

INDIR='tmp/' + CAN + '/intermediate_files/normalized/'
OUTDIR='tmp/' + CAN + '/intermediate_files/bpgs/'

INFILES = [f for f in listdir(INDIR) if isfile(join(INDIR, f))]
INFILES = ''.join(INFILES).replace('.P.rnk', ' ').split()

KEGG='data/c2.cp.kegg.v6.0.symbols.tsv'
GO='data/c5.bp.v6.0.symbols.tsv'

########################
### --- WORKFLOW --- ###
########################

### TRIGGER ###
rule all:
        input: expand(OUTDIR+'{infile}.bpg', infile=INFILES)

###  ###
rule p2bpg:
        input: INDIR+'{infile}.P.rnk'
        output: OUTDIR+'{infile}.bpg'
        run:
            shell("perl {PNL} {input} {GO} {KEGG} {PNL_MODE} > {output}")
