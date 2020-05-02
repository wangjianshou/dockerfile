#!/bin/bash
#第一个参数是物种名称，第二个参数是输出目录，第3个参数是参考基因组列表，空格隔开，最后放在双引号里，第四个参数是markers，可以使用通配符

strainphlan.py --ifn_samples $4 --ifn_markers $2/$1.markers.fasta --ifn_ref_genomes $3 --output_dir $2 --clades $1 --marker_in_clade 0.2
add_metadata_tree.py --ifn_trees $2/RAxML_bestTree.${1}.tree --ifn_metadatas $2/metadata.txt --metadatas group
plot_tree_graphlan.py --ifn_tree $2/RAxML_bestTree.${1}.tree.metadata --colorized_metadata group

