FROM continuumio/miniconda

MAINTAINER Lars Ridder <l.ridder@esciencecenter.nl> 
 
RUN conda config --set auto_update_conda false && \
conda install -y -c rdkit rdkit && \
conda install -y cython lxml nose coverage && \
conda clean -y -s -p -t -l -i

RUN pip install http://www.parallelpython.com/downloads/pp/pp-1.6.4.zip

RUN apt-get update && apt-get install  -y -q gcc 

ADD ./job /MAGMa/job

WORKDIR /MAGMa/job
RUN python setup.py develop

WORKDIR /root
RUN echo '[magma job]\nstructure_database.online = True\nstructure_database.service = http://www.emetabolomics.org/magma/molecules' > magma_job.ini

WORKDIR /data
VOLUME /data

ENTRYPOINT ["magma"]

