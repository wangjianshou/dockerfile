FROM r-base:3.5.0
RUN apt-get update && \
    apt-get install -y libcurl4-openssl-dev libxml2 libxml2-dev libssl-dev vim&& \
		echo "update.packages(checkBuilt=TRUE, ask=FALSE)" | R --slave && \
		echo "install.packages(\"BiocManager\");BiocManager::install(c(\"dada2\",\"DECIPHER\",\"phyloseq\"), version = \"3.8\")" | R --slave
    #echo "source(\"https://bioconductor.org/biocLite.R\");biocLite(\"dada2\")" | R --slave

RUN echo "PS1='${debian_chroot:+($debian_chroot)}\h:\w\n\$ '" >>~/.bashrc
COPY paired_dada2.R silva_nr_v132_train_set.fa.gz silva_species_assignment_v132.fa.gz SILVA_SSU_r132_March2018.RData /workdir/
