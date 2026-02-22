FROM ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive
ARG dependencies="wget curl unzip build-essential libz-dev \
                libncurses-dev libbz2-dev liblzma-dev git python3 gfortran \
                libreadline6-dev xorg-dev libblas-dev autotools-dev strace \
                gcc-multilib gobjc++ aptitude libpcre2-dev libcurl4 \
                libcurl4-openssl-dev default-jre default-jdk openjdk-8-jdk openjdk-8-jre \
                python3-pip bc pigz perl libdevel-size-perl \
                libfile-which-perl libtext-soundex-perl libjson-perl liburi-perl libwww-perl \
                gcc g++ libgomp1 python3-h5py cmake python3-xopen \
                cpio python3-setuptools bioperl libparallel-forkmanager-perl \
                zlib1g libpam-systemd tzdata"

RUN mkdir -p /Apps/ \
    && mkdir -p /opt \
    && apt-get -y update \
    && apt-get -y install ${dependencies} \
    && ln -s /usr/bin/python3 /usr/bin/python 

WORKDIR /Apps

#############################
###    Install RMBlast    ###
#############################

ENV RMBlast_Version=2.14.0

RUN mkdir RMBlast \
    && cd RMBlast \
    && wget -q ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${RMBlast_Version}/ncbi-blast-${RMBlast_Version}+-src.tar.gz \
    && wget -q  https://www.repeatmasker.org/rmblast/isb-${RMBlast_Version}+-rmblast.patch.gz \
    && tar xf ncbi-blast-${RMBlast_Version}+-src.tar.gz \
    && rm ncbi-blast-${RMBlast_Version}+-src.tar.gz \
    && gunzip isb-${RMBlast_Version}+-rmblast.patch.gz \
    && cd ncbi-blast-${RMBlast_Version}+-src \
    && patch -p1 < ../isb-${RMBlast_Version}+-rmblast.patch \
    && cd c++/ \
    && ./configure --with-mt --without-debug --without-krb5 --without-openssl --with-projects=scripts/projects/rmblastn/project.lst --prefix=/opt/rmblast \
    && make \
    && make install

ENV PATH="/opt/rmblast/bin:$PATH"

#############################
###    Install cd-hit   ###
#############################
ENV CDHIT_Version=4.8.1

RUN mkdir CDHIT \
    && cd CDHIT \
    && wget -q -O CDHIT.tar.gz https://github.com/weizhongli/cdhit/archive/refs/tags/V${CDHIT_Version}.tar.gz \
    && tar xf CDHIT.tar.gz \
    && rm CDHIT.tar.gz \
    && cd cdhit-${CDHIT_Version} \
    && make \
    && mkdir /opt/cd-hit \
    && PREFIX=/opt/cd-hit make install

ENV PATH="/opt/cd-hit:$PATH"

#############################
###     Install HMMER     ###
#############################
ENV HMMER_Version=3.4

RUN mkdir HMMER \
    && cd HMMER \
    && wget -q http://eddylab.org/software/hmmer/hmmer-${HMMER_Version}.tar.gz \
    && tar xf hmmer-${HMMER_Version}.tar.gz \
    && rm hmmer-${HMMER_Version}.tar.gz \
    && cd hmmer-${HMMER_Version} \
    && ./configure --prefix=/opt/hmmer \
    && make \
    && make install \
    && cd easel \
    && make install 

ENV PATH="/opt/hmmer/bin:$PATH"

#############################
###       Install R       ###
#############################
ENV R_Version=4.3.1

RUN curl -L https://cran.r-project.org/src/base/R-4/R-${R_Version}.tar.gz > R.tar.gz \
    && tar xf R.tar.gz \
    && cd R-${R_Version} \
    && ./configure --with-x=no \
    && make \
    && make check \
    && make install \
    && rm /Apps/R.tar.gz

#Install R-packages
ENV R_packages="'tidyverse'"
RUN R -e "install.packages( c( ${R_Packages} ), dependencies=TRUE, repos='http://cran.rstudio.com/')"

#############################
###    Install samtools   ###
#############################
ENV Samtools_Version=1.21

#Install htslib
RUN   wget https://github.com/samtools/htslib/releases/download/${Samtools_Version}/htslib-${Samtools_Version}.tar.bz2 \
      && tar xvjf htslib-${Samtools_Version}.tar.bz2 \
      && cd htslib-${Samtools_Version} \
      && ./configure \
      && make \
      && make install \
      && rm /Apps/htslib-${Samtools_Version}.tar.bz2

#Install samtools
RUN   wget https://github.com/samtools/samtools/releases/download/${Samtools_Version}/samtools-${Samtools_Version}.tar.bz2 \
      && tar -xvjf samtools-${Samtools_Version}.tar.bz2 \
      && cd samtools-${Samtools_Version} \
      && ./configure \
      && make \
      && make install \
      && rm /Apps/samtools-${Samtools_Version}.tar.bz2

#############################
###     Install MAFFT     ###
#############################
ENV MAFFT_Version=v7.526

RUN mkdir MAFFT \
    && cd MAFFT \
    && git clone https://gitlab.com/sysimm/mafft.git --branch ${MAFFT_Version} --single-branch \
    && cd ./mafft/core \
    && sed -i 's#^PREFIX =.*#PREFIX = /opt/mafft#' Makefile \
    && make \
    && make install

ENV PATH="/opt/mafft/bin:$PATH"

#############################
###    Install Bedtools   ###
#############################
ENV bedtools_Version=2.30.0

#Install Bedtools
RUN curl -L https://github.com/arq5x/bedtools2/archive/refs/tags/v${bedtools_Version}.tar.gz  > bedtools.tar.gz \
    && tar -zxvf bedtools.tar.gz \
    && cd bedtools2-${bedtools_Version}\
    && make \
    && make install \
    && rm /Apps/bedtools.tar.gz 

#############################
###      Install TRF      ###
#############################
ENV TRF_Version=4.09.1

RUN mkdir TRF \
    && cd TRF \
    && wget -q -O TRF.tar.gz https://github.com/Benson-Genomics-Lab/TRF/archive/refs/tags/v${TRF_Version}.tar.gz \
    && tar xf TRF.tar.gz \
    && rm TRF.tar.gz \
    && cd TRF-${TRF_Version} \
    && mkdir build \
    && cd build \
    && ../configure \
    && make \
    && mkdir /opt/trf \
    && cp ./src/trf /opt/trf/

ENV PATH="/opt/trf:$PATH"

#############################
###  Install RepeatScout  ###
#############################
ENV RepeatScout_Version=1.0.7

RUN mkdir RepeatScout \
    && cd RepeatScout \
    && wget -O RepeatScout.tar.gz -q https://github.com/Dfam-consortium/RepeatScout/archive/refs/tags/v${RepeatScout_Version}.tar.gz \
    && tar xf RepeatScout.tar.gz \
    && rm RepeatScout.tar.gz \
    && cd RepeatScout-${RepeatScout_Version} \
    && sed -i 's#^INSTDIR =.*#INSTDIR = /opt/RepeatScout#' Makefile \
    && sed -i 's#README#README.md#' Makefile \
    && make \
    && make install

ENV PATH="/opt/RepeatScout:$PATH"

#############################
###     Install RECON     ###
#############################
#&& wget -q https://www.repeatmasker.org/RepeatModeler/RECON-${RECON_Version}.tar.gz \
ENV RECON_Version=1.08

RUN mkdir RECON \
    && cd RECON \
    && wget -q -O RECON-${RECON_Version}.tar.gz https://blueprints.launchpad.net/ubuntu/+archive/primary/+sourcefiles/repeatmasker-recon/${RECON_Version}-8/repeatmasker-recon_${RECON_Version}.orig.tar.gz \
    && tar xf RECON-${RECON_Version}.tar.gz \
    && cd RECON-${RECON_Version} \
    && make -C src \
    && make -C src install \
    && cp 00README bin/ \
    && sed -i 's#^\$path =.*#$path = "/opt/RECON/bin";#' scripts/recon.pl \
    && cd /Apps/ \
    && mv /Apps/RECON/RECON-${RECON_Version} /opt/RECON

ENV PATH="/opt/RECON/bin:/opt/RECON/scripts:$PATH"

#############################
###  Install GenomeTools  ###
#############################
ENV GTools_Version=1.6.5

RUN mkdir genometools \
    && cd genometools \
    && wget -q -O gt.tar.gz https://github.com/genometools/genometools/archive/v${GTools_Version}.tar.gz \
    && tar xf gt.tar.gz \
    && rm gt.tar.gz \
    && cd genometools-${GTools_Version} \
    && make -j4 cairo=no \
    && make cairo=no prefix=/opt/genometools install

ENV PATH="/opt/genometools/bin:$PATH"

#############################
##    Install TESorter     ##
#############################
ENV TESorter_Version=1.4.7

RUN wget -O TESorter.tar.gz https://github.com/zhangrengang/TEsorter/archive/refs/tags/v${TESorter_Version}.tar.gz \
    && pip install TESorter.tar.gz \
    && rm TESorter.tar.gz

#############################
### Install LTR_retriever ###
#############################
ENV LTRRetriever_Version=3.0.1

RUN mkdir LTRRetriever \
    && cd LTRRetriever \
    && wget -q -O LTRRetriever.tar.gz https://github.com/oushujun/LTR_retriever/archive/refs/tags/v${LTRRetriever_Version}.tar.gz \
    && tar xf LTRRetriever.tar.gz \
    && rm LTRRetriever.tar.gz \
    && sed -i \
        -e 's#BLAST+=#BLAST+=/opt/rmblast/bin#' \
        -e 's#RepeatMasker=#RepeatMasker=/opt/RepeatMasker#' \
        -e 's#HMMER=#HMMER=/opt/hmmer/bin#' \
        -e 's#CDHIT=#CDHIT=/opt/cd-hit#' \
        LTR_retriever-${LTRRetriever_Version}/paths \
    && mv LTR_retriever-${LTRRetriever_Version} /opt/LTR_retriever

ENV PATH="/opt/LTR_retriever:$PATH"

#############################
###     Install NINJA     ###
#############################
ENV NINJA_Version=0.98

RUN mkdir NINJA \
    && cd NINJA \
    && wget -q -O NINJA.tar.gz https://github.com/TravisWheelerLab/NINJA/archive/refs/tags/${NINJA_Version}-cluster_only.tar.gz \
    && tar xf NINJA.tar.gz \
    && rm NINJA.tar.gz \
    && mv NINJA-${NINJA_Version}-cluster_only/* ./ \
    && rm -r NINJA-${NINJA_Version}-cluster_only/ \
    && cd NINJA \
    && make clean \
    && make all \
    && mv /Apps/NINJA /opt/NINJA

ENV PATH="/opt/NINJA/NINJA:$PATH"

#############################
###    Install UCSCdeps   ###
#############################
RUN mkdir /opt/ucsc_tools \
    && cd /opt/ucsc_tools \
    && wget -q https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit \
    && wget -q https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitInfo \
    && wget -q https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa \
    && chmod +x /opt/ucsc_tools/*

ENV PATH="/opt/ucsc_tools:$PATH"
COPY LICENSE.ucsc /opt/ucsc_tools/LICENSE

#############################
###     Install COSEG     ###
#############################
ENV COSEG_Version=0.2.4

RUN mkdir COSEG \
    && cd COSEG \
    && wget -O coseg-${COSEG_Version}.tar.gz -q https://github.com/rmhubley/coseg/archive/refs/tags/coseg-${COSEG_Version}.tar.gz \
    && tar xf coseg-${COSEG_Version}.tar.gz \
    && rm coseg-${COSEG_Version}.tar.gz \
    && cd coseg-coseg-${COSEG_Version} \
    && sed -i 's@#!.*perl@#!/usr/bin/perl@' preprocessAlignments.pl runcoseg.pl refineConsSeqs.pl \
    && sed -i 's#use lib "/usr/local/RepeatMasker";#use lib "/opt/RepeatMasker";#' preprocessAlignments.pl \
    && make \
    && mv /Apps/COSEG/coseg-coseg-${COSEG_Version} /opt/coseg

ENV PATH="/opt/coseg:$PATH"

#############################
###  Install RepeatMasker ###
#############################
ENV RepeatMasker_Version=4.1.8

RUN cd /opt \
    && wget -q https://www.repeatmasker.org/RepeatMasker/RepeatMasker-${RepeatMasker_Version}.tar.gz \
    && tar xf RepeatMasker-${RepeatMasker_Version}.tar.gz \
    && rm RepeatMasker-${RepeatMasker_Version}.tar.gz \
    && cd RepeatMasker \
    && perl configure \
       -hmmer_dir=/opt/hmmer/bin \
       -rmblast_dir=/opt/rmblast/bin \
       -libdir=/opt/RepeatMasker/Libraries \
       -trf_prgm=/opt/trf \
       -default_search_engine=rmblast

ENV PATH="/opt/RepeatMasker:/opt/RepeatMasker/util:$PATH"

#############################
### Install RepeatModeler ###
#############################
ENV RepeatModeler_Version=2.0.6a

ADD "https://www.random.org/cgi-bin/randbyte?nbytes=10&format=h" skipcache
RUN cd /opt \
    && wget -q -O RepeatModeler.tar.gz https://github.com/BioFalcon/RepeatModeler/archive/refs/tags/v${RepeatModeler_Version}.tar.gz \
    && tar xf RepeatModeler.tar.gz \
    && rm RepeatModeler.tar.gz \
    && mv RepeatModeler-${RepeatModeler_Version} RepeatModeler \
    && cd RepeatModeler \
    && perl configure \
         -cdhit_dir=/opt/cd-hit \
         -genometools_dir=/opt/genometools/bin \
         -ltr_retriever_dir=/opt/LTR_retriever \
         -mafft_dir=/opt/mafft/bin \
         -trf_dir=/opt/trf/ \
         -ninja_dir=/opt/NINJA/NINJA \
         -recon_dir=/opt/RECON/bin \
         -repeatmasker_dir=/opt/RepeatMasker \
         -rmblast_dir=/opt/rmblast/bin \
         -rscout_dir=/opt/RepeatScout \
         -ucsctools_dir=/opt/ucsc_tools

ENV PATH="/opt/RepeatModeler:/opt/RepeatModeler/util:$PATH"

#############################
###    Install seqkit     ###
#############################
ENV seqkit_Version=2.9.0

RUN mkdir /opt/seqkit \
    && wget -q -O seqkit.tar.gz https://github.com/shenwei356/seqkit/releases/download/v${seqkit_Version}/seqkit_linux_amd64.tar.gz \
    && tar xf seqkit.tar.gz \
    && rm seqkit.tar.gz \
    && mv seqkit /opt/seqkit

ENV PATH="/opt/seqkit/:${PATH}"

#############################
###  Install Python Deps  ###
#############################
ENV PythonLibs="tqdm==4.65.0 pandas==2.0.3  portion==2.4.1 \
                matplotlib==3.7.1 numpy==1.26.4 scikit-image==0.21.0 \
                biopython==1.81 configparser==6.0.1 dendropy==4.5.2 \
                psutil==5.9.8 scipy==1.12.0 pybedtools==0.9.1 \
                pysam==0.22.0 polars==1.5.0 seaborn==0.13.2 \
                xopen==2.0.2 drmaa==0.7.6"

RUN pip3 install ${PythonLibs}

#############################
##     Install bedops     ###
#############################
ENV BEDOPS_Version=2.4.41

RUN wget -O bedops.tar.gz https://github.com/bedops/bedops/archive/v${BEDOPS_Version}.tar.gz \
    && tar xf bedops.tar.gz \ 
    && cd bedops-2.4.41 \
    && make \
    && make install \
    && mkdir /opt/bedops \
    && cp bin/* /opt/bedops
    
ENV PATH="/opt/bedops/:${PATH}"

#############################
###         Cleanup       ###
#############################
RUN rm -r /Apps 
