FROM --platform=linux/amd64 rocker/verse

## Install only the R packages we need for the project
RUN install2.r --error \
    Seurat \
    patchwork
