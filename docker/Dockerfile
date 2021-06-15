FROM debian:buster-20210511-slim

MAINTAINER Matt Eldridge "matthew.eldridge@cruk.cam.ac.uk"

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update
RUN apt-get install -y apt-utils

# update locale
RUN apt-get install -y locales
RUN sed -i 's/^# en_GB/en_GB/' /etc/locale.gen
RUN locale-gen

# install packages and libraries needed to build/install Shiny server, R and
# R packages/dependencies
RUN apt-get install -y build-essential
RUN apt-get install -y libxml2-dev
RUN apt-get install -y libssl-dev
RUN apt-get install -y libreadline7
RUN apt-get install -y libcurl4-openssl-dev
RUN apt-get install -y curl
RUN apt-get install -y gdebi-core
RUN apt-get clean

# install Shiny server
ARG SHINY_SERVER_VERSION=1.5.16.958-amd64
RUN curl -O https://download3.rstudio.org/ubuntu-14.04/x86_64/shiny-server-${SHINY_SERVER_VERSION}.deb
RUN gdebi -n shiny-server-${SHINY_SERVER_VERSION}.deb
RUN rm shiny-server-${SHINY_SERVER_VERSION}.deb

# expose port that Shiny server listens on (can be mapped to another port on the host)
EXPOSE 3838

# ensure the user account that is running the shiny server process has write privilege for /var/lib/shiny-server
RUN chmod ugo+rwx /var/lib/shiny-server

# command to run when starting the container using 'docker run'
CMD exec shiny-server >> /var/log/shiny-server/shiny-server.log 2>&1

# install R
RUN apt-key adv --keyserver keys.gnupg.net --recv-key E19F5F87128899B192B1A2C2AD5F960A256A04AF
RUN echo 'deb http://cloud.r-project.org/bin/linux/debian buster-cran40/' >> /etc/apt/sources.list
RUN apt-get update
RUN apt-get install -y -t buster-cran40 r-base

# update pre-installed R packages
RUN R -e 'update.packages(ask = FALSE)'

# install R packages
RUN R -e 'install.packages("tidyverse")'
RUN R -e 'install.packages("optparse")'
RUN R -e 'install.packages("shiny")'
RUN R -e 'install.packages("shinyjs")'
RUN R -e 'install.packages("colourpicker")'
RUN R -e 'install.packages("DT")'

# install Bioconductor packages
RUN R -e 'install.packages("BiocManager")'
RUN R -e 'BiocManager::install("Biobase")'
RUN R -e 'BiocManager::install("QDNAseq")'

# install rascal package
ARG RASCAL_VERSION=0.6.0
COPY rascal_${RASCAL_VERSION}.tgz /opt
RUN R CMD INSTALL /opt/rascal_${RASCAL_VERSION}.tgz

# extract rascal Shiny app and create link in shiny server deployment directory
RUN tar zxf /opt/rascal_${RASCAL_VERSION}.tgz -C /opt
RUN ln -s /opt/rascal/shiny /srv/shiny-server/rascal

