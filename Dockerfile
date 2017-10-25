## Start with the tidyverse docker image
FROM rocker/tidyverse:latest

MAINTAINER "Sam Abbott" sam.abbott@bristol.ac.uk

ADD . /home/seabbs

RUN Rscript -e 'devtools::install_dev_deps("/home/seabbs")'

RUN Rscript -e 'devtools::install_github("hadley/pkgdown")'
