FROM continuumio/miniconda

MAINTAINER Lars Ridder <l.ridder@esciencecenter.nl> 
 
RUN /opt/conda/bin/conda install -y -q -c https://conda.anaconda.org/rdkit rdkit && \
/opt/conda/bin/conda install -y cython lxml nose coverage && \
/opt/conda/bin/conda clean -y -s -p -t -l -i

ENV PATH /opt/conda/bin:$PATH 

RUN pip install http://www.parallelpython.com/downloads/pp/pp-1.6.4.zip

RUN apt-get update && apt-get install  -y -q gcc 

RUN git clone -b master https://github.com/NLeSC/MAGMa.git

WORKDIR /MAGMa/job
RUN /opt/conda/bin/python setup.py develop

WORKDIR /
RUN echo '[magma job]\nstructure_database.online = True\nstructure_database.service = http://www.emetabolomics.org/magma/molecules' > magma_job.ini

WORKDIR /data
VOLUME /data

ENTRYPOINT ["magma"]
