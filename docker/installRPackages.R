#20161110 use following command to install R packages to all systems on computing nodes
#root@n100:~/install# ~/bin/runOnAllNodes.sh /usr/bin/Rscript ~/installRPackages.R
pkgNames = c("ggplot2", "gridExtra", "magrittr", "gplots");
for (n in pkgNames){
	install.packages(n, repos="http://cran.us.r-project.org")
	#on simm cluster
	#install.packages(n, repos="http://cran.us.r-project.org", lib="/y/program/lib/R/site-library")
}
