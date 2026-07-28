ARG R_VERSION=4.5.1
FROM rocker/shiny-verse:${R_VERSION}

LABEL org.opencontainers.image.title="SimpleMicrobiome"
LABEL org.opencontainers.image.description="R Shiny application for microbiome analysis"
LABEL org.opencontainers.image.source="https://github.com/yjcho2252/SimpleMicrobiome"

ENV DEBIAN_FRONTEND=noninteractive
ENV R_DEFAULT_INTERNET_TIMEOUT=600
ENV SIMPLEMICROBIOME_R_LIBRARY=/usr/local/lib/R/site-library

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

# Restore exactly the package versions and GitHub commits recorded from the
# tested image. install.R also verifies the critical runtime packages.
RUN Rscript install.R \
&& cd / \
&& rm -rf /tmp/simplemicrobiome-dependencies

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
