#Make a base image without accucopy

all:
	date
	ls pyflow
	ls bashrc
	ls samtools
	ls strelka-2.9.10.centos6_x86_64
	docker login
	docker pull ubuntu:18.04
	docker build -t polyactis/accucopy:18 .
	docker push polyactis/accucopy:18
	date
