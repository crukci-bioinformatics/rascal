#!/bin/bash

sed "s/run_as USER/run_as $USER/" shiny-server.conf.template > shiny-server.conf

dir=`dirname $0`

mkdir -p ${dir}/logs
mkdir -p ${dir}/lib

singularity run \
	--cleanenv \
	--bind ${dir}/shiny-server.conf:/etc/shiny-server/shiny-server.conf \
	--bind ${dir}/lib:/var/lib/shiny-server \
	--bind ${dir}/logs:/var/log/shiny-server \
	${dir}/rascal.sif

# add this line to deploy local changes to the Shiny app during development
#	--bind ${dir}/../inst/shiny:/srv/shiny-server/rascal \

