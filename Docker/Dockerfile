FROM ghcr.io/prefix-dev/pixi:0.78.0

ARG IMGT_HLA_VERSION=""
ARG ARCASHLA_REF_DIR=/opt/arcashla-ref
ENV ARCASHLA_REF_DIR=${ARCASHLA_REF_DIR}


WORKDIR /app/arcasHLA
COPY . /app/arcasHLA

# Installs git, git-lfs, pigz and the rest of the runtime dependencies
# from pixi.toml/pixi.lock into the "default" environment.
RUN pixi install --frozen --environment default

# Create the shell-hook bash script to activate the environment
RUN pixi shell-hook > /shell-hook.sh

# Ensure the arcasHLA main script gets found. This must be appended to the
# shell-hook script (not set via ENV) because the pixi shell-hook hardcodes
# an absolute PATH value that would otherwise overwrite/clobber an ENV PATH.
RUN echo 'export PATH="/app/arcasHLA:$PATH"' >> /shell-hook.sh

# extend the shell-hook script to run the command passed to the container
RUN echo 'exec "$@"' >> /shell-hook.sh

RUN if [ -n "${IMGT_HLA_VERSION}" ]; then \
        git clone --depth 1 --branch "${IMGT_HLA_VERSION}" \
            https://github.com/ANHIG/IMGTHLA.git /opt/imgt-hla && \
        ./arcasHLA reference build \
            --imgt /opt/imgt-hla \
            --outdir "${ARCASHLA_REF_DIR}" \
            --version "${IMGT_HLA_VERSION}"; \
    fi

# set the entrypoint to the shell-hook script (activate the environment and run the command)
# no more pixi needed in the prod container
ENTRYPOINT ["/bin/bash", "/shell-hook.sh"]

CMD ["arcasHLA"]
