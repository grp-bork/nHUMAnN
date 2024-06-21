FROM ubuntu:22.04

LABEL maintainer="cschu1981@gmail.com"
LABEL version="0.7"
LABEL description="This is a Docker Image for the reCOGnise tool."


ARG DEBIAN_FRONTEND=noninteractive

RUN apt update
RUN apt upgrade -y

RUN apt-get install -y wget python3-pip git dirmngr gnupg ca-certificates build-essential libssl-dev libcurl4-gnutls-dev libxml2-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev rsync bowtie2 diamond-aligner

RUN apt clean



RUN wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
  bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/software/miniconda3 && \
  rm -f Miniconda3-latest-Linux-x86_64.sh 


ENV PATH /opt/software/miniconda3/bin:$PATH

RUN conda install -c conda-forge -c bioconda -y 'diamond==2.0.15'

RUN git clone https://github.com/biobakery/MetaPhlAn.git metaphlan && \
  cd metaphlan && \
  git checkout 4.0.6 && \
  pip install .

RUN git clone https://github.com/biobakery/humann.git && \
  cd humann && \
  git checkout e3b05c5e06220c991cb7f7578315ca495492224d && \
  sed -i "363 s/3/4/" humann/config.py && \
  sed -i "205 s/\$/.strip()/" humann/search/prescreen.py && \
  sed -i '153i\                print(version_found)' humann/search/prescreen.py && \
  sed -i "100 s/message/line/" humann/search/prescreen.py && \
  sed -i "727 s/not config.metaphlan_v3_db_matching_uniref in file/False/" humann/humann.py && \
  sed -i "796 s/not config.matching_uniref in file/False/" humann/humann.py && \
  python setup.py install
