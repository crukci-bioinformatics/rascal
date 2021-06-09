FROM debian:buster-20210511-slim

MAINTAINER Matt Eldridge "matthew.eldridge@cruk.cam.ac.uk"

RUN apt-get update
RUN apt-get install -y build-essential
RUN apt-get install -y libxml2-dev
RUN apt-get install -y libssl-dev
RUN apt-get install -y libreadline7
RUN apt-get install -y curl
RUN apt-get install -y gdebi-core
RUN apt-get clean

# install R
ARG R_VERSION=4.1.0
RUN curl -O https://cdn.rstudio.com/r/debian-10/pkgs/r-${R_VERSION}_1_amd64.deb
RUN gdebi -n r-${R_VERSION}_1_amd64.deb
RUN rm r-${R_VERSION}_1_amd64.deb
RUN ln -s /opt/R/${R_VERSION}/bin/R /usr/local/bin/R
RUN ln -s /opt/R/${R_VERSION}/bin/Rscript /usr/local/bin/Rscript

# install Shiny server
ARG SHINY_SERVER_VERSION=1.5.16.958-amd64
RUN curl -O https://download3.rstudio.org/ubuntu-14.04/x86_64/shiny-server-${SHINY_SERVER_VERSION}.deb
RUN gdebi -n shiny-server-${SHINY_SERVER_VERSION}.deb
RUN rm shiny-server-${SHINY_SERVER_VERSION}.deb

# install R packages
RUN R -e 'install.packages("shiny", repos = "https://cloud.r-project.org")'
RUN R -e 'install.packages("tidyverse", repos = "https://cloud.r-project.org")'
RUN R -e 'install.packages("DT", repos = "https://cloud.r-project.org")'
RUN R -e 'install.packages("shinyjs", repos = "https://cloud.r-project.org")'
RUN R -e 'install.packages("colourpicker", repos = "https://cloud.r-project.org")'

# install Bioconductor packages
RUN R -e 'install.packages("BiocManager", repos = "https://cloud.r-project.org")'
RUN R -e 'BiocManager::install("Biobase")'
RUN R -e 'BiocManager::install("QDNAseq")'

# install rascal package
ARG RASCAL_VERSION=0.6.0
COPY rascal_${RASCAL_VERSION}.tgz /opt
RUN R CMD INSTALL /opt/rascal_${RASCAL_VERSION}.tgz

# extract rascal Shiny app and create link in shiny server deployment directory
RUN tar zxf /opt/rascal_${RASCAL_VERSION}.tgz -C /opt
RUN ln -s /opt/rascal/shiny /srv/shiny-server/rascal

# ensure the user account that is running the shiny server process has write privilege for /var/lib/shiny-server
RUN chmod ugo+rwx /var/lib/shiny-server

# expose port that Shiny server listens on (can be mapped to another port on the host)
EXPOSE 3838

# command to run when starting the container using 'docker run'
CMD exec shiny-server >> /var/log/shiny-server/shiny-server.log 2>&1

