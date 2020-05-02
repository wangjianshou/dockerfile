#!/bin/bash
run_metaphlan2.py -i $1/runid -o $1/step1.out/ -p 10 --strainphlan
strainphlan.py --ifn_samples $1/step1.out/*markers --output_dir $1/step2.out --print_clades_only >$1/step2.out/clades.txt
for i in `cat $1/step2.out/clades.txt`
do
  strainphlan.py --ifn_samples $1/step1.out/*markers --output_dir $1/step2.out --clades $i --marker_in_clade 0.2 --nprocs_main 4
done

<<comment
for i in `find $1/step2.out -name RAxML_bestTree.*.tree`
do
  add_metadata_tree.py --ifn_trees $i --ifn_metadatas $1/metadata.txt --metadatas group
done
comment
