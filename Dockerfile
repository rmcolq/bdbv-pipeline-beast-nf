FROM condaforge/mambaforge:latest AS conda

COPY . ./bdbv-pipeline-beast/

RUN /opt/conda/bin/mamba env create -f ./bdbv-pipeline-beast/environment.yml

ENV PATH=/opt/conda/envs/bdbv-pipeline-beast/bin:$PATH

ENV MPLCONFIGDIR="."

CMD ["/bin/bash"]
