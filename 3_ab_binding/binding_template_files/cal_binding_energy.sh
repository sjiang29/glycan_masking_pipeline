#!/bin/bash
 
time /sb/meilerapps/rosetta/rosetta-3.13/main/source/bin/rosetta_scripts.default.linuxgccrelease -l $1 -parser:protocol template_docking_analysis.xml -out:file:score_only -out:file:scorefile InAn_cleaned_$2.sc > InAn.log 
