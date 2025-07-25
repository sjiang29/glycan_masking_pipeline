#!/bin/bash
 
time /sb/meilerapps/rosetta/rosetta-3.13/main/source/bin/rosetta_scripts.default.linuxgccrelease -l $1 -parser:protocol $2 @AddChain_FastRelax_Multi.options -nstruct 1 -out:prefix Add_Rlx_$3 > AddRel.log
