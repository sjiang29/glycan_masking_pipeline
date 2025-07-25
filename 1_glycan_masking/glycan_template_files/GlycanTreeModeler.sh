#!/bin/bash
time /dors/meilerlab/apps/rosetta/rosetta-3.13/main/source/bin/rosetta_scripts.default.linuxgccrelease -s $1 -parser:protocol GlycanTreeModeler.xml -nstruct $2 @GlycanTreeModeler.options -out:prefix Glyc_ -out:file:scorefile Glyc_score.sc > Glyc.log
