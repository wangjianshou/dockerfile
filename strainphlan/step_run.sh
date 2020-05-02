docker run --rm -v /k11e/:/k11e/ wangjianshou/strainphlan:v2 run_metaphlan2.py -i `pwd`/runid -o `pwd`/step1.out/ -p 10 --strainphlan

docker run --rm -v `pwd`:/work wangjianshou/strainphlan:v2 bash -c "strainphlan.py --ifn_samples /work/step1.out/*markers --output_dir /work/step2.out --print_clades_only >/work/step2.out/clades.txt"

for i in `cat step2.out/clades.txt`
do
  docker run --rm -v `pwd`:/work wangjianshou/strainphlan:v2 strainphlan.py --ifn_samples /work/step1.out/*markers --output_dir /work/step2.out --clades $i --marker_in_clade 0.2 --nprocs_main 10
done

<<comment
for i in `find step2.out -name RAxML_bestTree.*.tree`
do
  docker run --rm -v `pwd`:/work wangjianshou/strainphlan:v2 add_metadata_tree.py --ifn_trees $i --ifn_metadatas metadata.txt --metadatas group
done
comment
