#!/bin/bash

mkdir -p logs
chmod ugo+rwx logs

docker run \
	--user shiny \
	--name rascal \
	--rm \
	--detach \
	--publish 8080:3838 \
	--volume ${PWD}/shiny-server.conf:/etc/shiny-server/shiny-server.conf \
	--volume ${PWD}/logs:/var/log/shiny-server \
	crukcibioinformatics/rascal

# add this line to deploy local changes to the Shiny app during development
#	--volume ${PWD}/../inst/shiny:/srv/shiny-server/rascal \

