
version ?= latest

build:
	docker build --tag crukcibioinformatics/rascal:${version} .

release: build
	docker push crukcibioinformatics/rascal:${version}

