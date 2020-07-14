# Add Accucopy to the base image

tag_name = allgenomes

all:
	date
	docker login
	docker login docker.yfish.x
	-mkdir binary
	cd ../../src_o/; make debug;
	rm -rf ./binary/*
	cp -apr ../../src_o/debug/* ./binary/
	cp -apr configure ./binary/
	docker pull polyactis/accucopy:18
	docker build -t polyactis/accucopy:latest .
	docker tag polyactis/accucopy:latest polyactis/accucopy:$(tag_name)
	#docker commit -m "purity (-0.1,0) be set to 0.03, (1,1.1) set to 0.99" -a "Yu Huang" d09dd0f24b58-Container-ID polyactis/accucopy:max_purity_1
	docker push polyactis/accucopy:$(tag_name)
	docker push polyactis/accucopy:latest
	docker tag polyactis/accucopy:latest docker.yfish.x/accucopy:$(tag_name)
	docker tag docker.yfish.x/accucopy:$(tag_name) docker.yfish.x/accucopy:latest
	docker push docker.yfish.x/accucopy
	date
