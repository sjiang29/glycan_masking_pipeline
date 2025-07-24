#!/bin/bash
time /dors/meilerlab/apps/rosetta/rosetta-3.13/main/source/bin/rosetta_scripts.default.linuxgccrelease -s $1 -parser:protocol FastDesign.xml @FastDesign.options -nstruct 1 -out:file:scorefile FD_score.sc -out:prefix Des_> FastDesign.log


