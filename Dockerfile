FROM ubuntu:26.04

LABEL maintainer="cschu1981@gmail.com"
LABEL version="3.9.cs.1x"
LABEL description="This is a Docker image for patching HUMAnN 3.9."


ARG DEBIAN_FRONTEND=noninteractive

RUN apt update
RUN apt upgrade -y

RUN apt-get install -y wget python3-pip python-is-python3 git dirmngr gnupg ca-certificates build-essential libssl-dev libcurl4-gnutls-dev libxml2-dev rsync
#bowtie2 diamond-aligner

RUN apt clean

RUN mkdir -p /opt/software
WORKDIR /opt/software

RUN git clone https://github.com/biobakery/MetaPhlAn.git metaphlan && \
  cd metaphlan && \
  git checkout 4.0.6 && \
  pip install --break-system-packages .

RUN git clone https://github.com/biobakery/humann.git && \
  cd humann && \
  git checkout v3.9 && \
  sed -i "s/\.readline()/.readline().rstrip()/g" humann/search/prescreen.py && \
  # git checkout e3b05c5e06220c991cb7f7578315ca495492224d && \
  # sed -i "363 s/3/4/" humann/config.py && \
  # sed -i "205 s/\$/.strip()/" humann/search/prescreen.py && \
  # sed -i '153i\                print(version_found)' humann/search/prescreen.py && \
  # sed -i "100 s/message/line/" humann/search/prescreen.py && \
  # sed -i "727 s/not config.metaphlan_v3_db_matching_uniref in file/False/" humann/humann.py && \
  # sed -i "796 s/not config.matching_uniref in file/False/" humann/humann.py && \
  # sed -i "s/vOct22/vJun23/" humann/config.py && \
  #python setup.py install
  pip install --break-system-packages .



#RUN git clone https://github.com/biobakery/humann.git && \
#  cd humann && \
#  git checkout e3b05c5e06220c991cb7f7578315ca495492224d && \
#  sed -i "363 s/3/4/" humann/config.py && \
#  sed -i "205 s/\$/.strip()/" humann/search/prescreen.py && \
#  sed -i '153i\                print(version_found)' humann/search/prescreen.py && \
#  sed -i "100 s/message/line/" humann/search/prescreen.py && \
#  sed -i "727 s/not config.metaphlan_v3_db_matching_uniref in file/False/" humann/humann.py && \
#  sed -i "796 s/not config.matching_uniref in file/False/" humann/humann.py && \
#  python setup.py install
