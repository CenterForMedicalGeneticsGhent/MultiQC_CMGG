ARG MULTIQC_VERSION=v1.33

FROM ghcr.io/multiqc/multiqc:${MULTIQC_VERSION}

LABEL org.opencontainers.image.title="MultiQC_CMGG" \
    org.opencontainers.image.description="MultiQC plugin for CMGG" \
    org.opencontainers.image.source="https://github.com/CenterForMedicalGeneticsGhent/MultiQC_CMGG" \
    org.opencontainers.image.licenses="MIT"

COPY . /src/

# Install dependencies and the package
RUN pip install --no-cache-dir /src/
