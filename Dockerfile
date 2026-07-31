# syntax=docker/dockerfile:1.7
ARG R_VERSION=4.5.1
FROM rocker/shiny-verse:${R_VERSION}

LABEL org.opencontainers.image.title="SimpleMicrobiome"
LABEL org.opencontainers.image.description="R Shiny application for microbiome analysis"
LABEL org.opencontainers.image.source="https://github.com/yjcho2252/SimpleMicrobiome"

ENV DEBIAN_FRONTEND=noninteractive
ENV R_DEFAULT_INTERNET_TIMEOUT=600
ENV SIMPLEMICROBIOME_R_LIBRARY=/usr/local/lib/R/site-library
ENV RENV_PATHS_CACHE=/renv/cache
ENV RENV_CONFIG_CACHE_SYMLINKS=FALSE

RUN apt-get update && apt-get install -y --no-install-recommends \
build-essential \
ca-certificates \
cmake \
curl \
gfortran \
git \
libblas-dev \
libcurl4-openssl-dev \
libfontconfig1-dev \
libfreetype6-dev \
libfftw3-dev \
libfribidi-dev \
libgdal-dev \
libgeos-dev \
libgit2-dev \
libglpk-dev \
libgmp3-dev \
libgsl-dev \
libharfbuzz-dev \
libjpeg-dev \
liblapack-dev \
libmpfr-dev \
libnlopt-dev \
libpng-dev \
libproj-dev \
libssl-dev \
libtiff-dev \
libudunits2-dev \
libxml2-dev \
libxt-dev \
make \
pandoc \
pkg-config \
zlib1g-dev \
&& rm -rf /var/lib/apt/lists/*

WORKDIR /tmp/simplemicrobiome-dependencies

COPY install.R renv.lock ./

# Reuse downloaded/built packages across Docker builds and retry transient
# repository failures without discarding packages installed earlier in the run.
RUN --mount=type=cache,id=simplemicrobiome-renv-r451,target=/renv/cache,sharing=locked \
    set -eu; \
    attempt=1; \
    while true; do \
      if Rscript --vanilla install.R; then \
        break; \
      fi; \
      if [ "$attempt" -ge 5 ]; then \
        echo "install.R failed after ${attempt} attempts"; \
        exit 1; \
      fi; \
      delay=$((attempt * 15)); \
      echo "install.R failed; retrying in ${delay} seconds"; \
      sleep "$delay"; \
      attempt=$((attempt + 1)); \
    done

RUN cd / && rm -rf /tmp/simplemicrobiome-dependencies

RUN rm -rf /srv/shiny-server/*

WORKDIR /srv/shiny-server

COPY --chown=shiny:shiny app.R ./app.R
COPY --chown=shiny:shiny modules ./modules
COPY --chown=shiny:shiny sample ./sample
COPY --chown=shiny:shiny www ./www
COPY --chown=shiny:shiny docs ./docs

EXPOSE 3838

HEALTHCHECK --interval=30s --timeout=5s --start-period=60s --retries=3 \
CMD curl --fail http://localhost:3838/ || exit 1

CMD ["/usr/bin/shiny-server"]
