FROM condaforge/miniforge3:24.7.1-0

LABEL maintainer="cschu1981@gmail.com"
LABEL version="2.2"
LABEL description="This is a Docker image for humann3"

ARG DEBIAN_FRONTEND=noninteractive

RUN apt update
RUN apt upgrade -y

RUN apt install -y wget python3-pip git
RUN apt clean

WORKDIR /opt/software
RUN	git clone https://github.com/biobakery/MetaPhlAn.git metaphlan && \
  	cd metaphlan && \
  	git checkout 4.0.6 && \
	python -m pip install . && \
	cd - && \
	rm -rf metaphlan && \
	git clone https://github.com/biobakery/humann.git

RUN conda create -y -n mp4_vOct22_v37 \
	"python<3.12" \
	# && cd /opt/software/metaphlan \
	# && /opt/conda/envs/mp4_vOct22_v37/bin/python -m pip install . \
	&& cd /opt/software/humann \ 
	&& git checkout e3b05c5e06220c991cb7f7578315ca495492224d \
	&& sed -i "363 s/3/4/" humann/config.py \
	&& sed -i "205 s/\$/.strip()/" humann/search/prescreen.py \
	&& sed -i '153i\                print(version_found)' humann/search/prescreen.py \
	&& sed -i "100 s/message/line/" humann/search/prescreen.py \
	&& sed -i "727 s/not config.metaphlan_v3_db_matching_uniref in file/False/" humann/humann.py \
	&& sed -i "796 s/not config.matching_uniref in file/False/" humann/humann.py \
	&& /opt/conda/envs/mp4_vOct22_v37/bin/python setup.py install \
	&& git reset --hard \
	&& conda clean -afy \
	&& printf '#!/bin/sh\nexec /opt/conda/envs/mp4_vOct22_v37/bin/humann "$@"\n' > /usr/local/bin/humann_37 \
	&& printf '#!/bin/sh\nexec /opt/conda/envs/mp4_vOct22_v37/bin/humann_reduce_table "$@"\n' > /usr/local/bin/humann_reduce_table_37 \
	&& printf '#!/bin/sh\nexec /opt/conda/envs/mp4_vOct22_v37/bin/humann_renorm_table "$@"\n' > /usr/local/bin/humann_renorm_table_37 \
	&& printf '#!/bin/sh\nexec /opt/conda/envs/mp4_vOct22_v37/bin/humann_split_stratified_table "$@"\n' > /usr/local/bin/humann_split_stratified_table_37 \
	&& printf '#!/bin/sh\nexec /opt/conda/envs/mp4_vOct22_v37/bin/humann_join_tables "$@"\n' > /usr/local/bin/humann_join_tables_37 \
	&& chmod +x /usr/local/bin/humann*_37


RUN conda create -y -n mp4_vJun23_v39 \
 	"python<3.12" \
	# && cd /opt/software/metaphlan \
	# && /opt/conda/envs/mp4_vJun23_v39/bin/python -m pip install . \
	&& cd /opt/software/humann \
	&& git checkout v3.9 \
	&& sed -i "s/\.readline()/.readline().rstrip()/g" humann/search/prescreen.py \
	&& /opt/conda/envs/mp4_vJun23_v39/bin/python setup.py install \
	&& cd - && rm -rf /opt/software/humann \
 	&& conda clean -afy \
	&& printf '#!/bin/sh\nexec /opt/conda/envs/mp4_vJun23_v39/bin/humann "$@"\n' > /usr/local/bin/humann_39 \
	&& printf '#!/bin/sh\nexec /opt/conda/envs/mp4_vJun23_v39/bin/humann_reduce_table "$@"\n' > /usr/local/bin/humann_reduce_table_39 \
	&& printf '#!/bin/sh\nexec /opt/conda/envs/mp4_vJun23_v39/bin/humann_renorm_table "$@"\n' > /usr/local/bin/humann_renorm_table_39 \
	&& printf '#!/bin/sh\nexec /opt/conda/envs/mp4_vJun23_v39/bin/humann_split_stratified_table "$@"\n' > /usr/local/bin/humann_split_stratified_table_39 \
	&& printf '#!/bin/sh\nexec /opt/conda/envs/mp4_vJun23_v39/bin/humann_join_tables "$@"\n' > /usr/local/bin/humann_join_tables_39 \
	&& chmod +x /usr/local/bin/humann*_39
	


#    && printf '#!/bin/sh\nexec /opt/conda/envs/genomespot/bin/python "$@"\n' \
#        > /usr/local/bin/python_gs \ 
#    && chmod +x /usr/local/bin/python_gs
