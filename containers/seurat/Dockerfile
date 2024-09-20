# Use the official R image with version 4.4.1
FROM r-base:4.4.1

# Install system dependencies for R packages
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    zlib1g-dev \
    libhdf5-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e "install.packages('remotes')"
RUN R -e "install.packages(c('argparse', 'dplyr', 'ggplot2', 'seurat'))"