#!/bin/bash
#第一个参数是样本第一步输出的markers文件，使用通配符指定所有样本的markers文件，第二个参数是输出目录
db="/home/linuxbrew/.linuxbrew/Cellar/strainphlan/2.6.0"
strainphlan.py --ifn_samples $1 --output_dir $2 --print_clades_only > $2/clades.txt
for i in `cat ${2}/clades.txt`
do
  [ -d "$2/$i" ] || mkdir -p "$2/$i"
	extract_markers.py --mpa_pkl ${db}/db_v20/mpa_v20_m200.pkl --ifn_markers ${db}/db_markers/all_markers.fasta --clade "$i" --ofn_markers "$2/$i/$i.markers.fasta" >/dev/null
  grep "$i" ${db}/species2genomes.txt >"$2/$i/$i.Genome.txt"
done

#strainphlan.py --ifn_samples *.markers --ifn_markers ${db}/db_markers/${i}.markers.fasta  --output_dir $2/$i --clades $i --marker_in_clade 0.2

