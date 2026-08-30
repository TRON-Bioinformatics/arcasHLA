FROM continuumio/miniconda3:latest

WORKDIR /app

COPY environment.yml /tmp/environment.yml
RUN conda env create --file /tmp/environment.yml && \
    conda clean --all --yes

COPY . .

ENV PATH="/app:/opt/conda/envs/arcas-hla/bin:${PATH}"
