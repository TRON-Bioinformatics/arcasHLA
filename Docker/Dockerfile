FROM continuumio/miniconda3:latest

WORKDIR /app

COPY environment.yml /tmp/environment.yml
RUN conda env create --file /tmp/environment.yml && \
    conda clean --all --yes

COPY . .

ENV PATH="/opt/conda/envs/arcas-hla/bin:/app:${PATH}"
