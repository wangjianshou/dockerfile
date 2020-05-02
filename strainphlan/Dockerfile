FROM biobakery/strainphlan:2.7.7
COPY mpa_v20_m200.fna mpa_v20_m200.pkl /opt/conda/bin/metaphlan_databases/
RUN bowtie2-build /opt/conda/bin/metaphlan_databases/mpa_v20_m200.fna /opt/conda/bin/metaphlan_databases/mpa_v20_m200 
RUN echo "PS1='${debian_chroot:+($debian_chroot)}\h:\w\n\$ '" >>~/.bashrc
RUN apt-get update && apt-get -y install vim
WORKDIR /work
COPY run_metaphlan2.py run.sh /usr/bin/
#RUN hg clone https://hg@bitbucket.org/nsegata/graphlan
#RUN ln -s /home/linuxbrew/graphlan/graphlan.py /home/linuxbrew/.linuxbrew/bin/ && ln -s /home/linuxbrew/graphlan/graphlan_annotate.py /home/linuxbrew/.linuxbrew/bin/
#RUN chmod +x /home/linuxbrew/step1_strainphlan.sh /home/linuxbrew/step2_strainphlan.sh
#RUN apt-get update && apt-get -y install wget
#RUN wget -O /usr/local/bin/gosu "https://github.com/tianon/gosu/releases/download/1.7/gosu-amd64" \
#    && chmod +x /usr/local/bin/gosu \
#		&& apt-get -y remove wget \
#		&& gosu nobody true
