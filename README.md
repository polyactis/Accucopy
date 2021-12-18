# Accucopy

Accucopy is a CNA-calling method that extends our previous Accurity model to predict both total (TCN) and allele-specific copy numbers (ASCN) for the tumor genome. Accucopy adopts a tiered Gaussian mixture model coupled with an innovative autocorrelation-guided EM algorithm to find the optimal solution quickly. The Accucopy model utilizes information from both total sequencing coverage and allelic sequencing coverage. Through comparative analyses in both simulation and real-sequencing samples, we demonstrate that Accucopy is more accurate than existing methods

## Publication

X Fan, G Luo, YS Huang# (2021) BMC Bioinformatics. Accucopy: Accurate and Fast Inference of Allele-specific Copy Number Alterations from Low-coverage Low-purity Tumor Sequencing Data.


# Compile

To compile, go into src_o/ and type "make debug" or "make release".

WARNING：Accucopy is not easy to compile. Check https://www.yfish.org/display/PUB/Accucopy#Accucopy-3.5Compilesourcecode(foradvancedusers) for how to compile.

If you just want to use Accucopy, using its docker version would be much easier: https://www.yfish.org/display/PUB/Accucopy#Accucopy-Docker.

You can also download the pre-compiled release. Check https://www.yfish.org/display/PUB/Accucopy for dependent packages to run Accucopy.

# Reference genome
We provide two different versions of human reference genomes, hs37d5 and hs38d1, downloadable from https://www.yfish.org/display/PUB/Accucopy#Accucopy-3.6Downloadareferencegenomefolder.

We recommend users to re-align reads against one of our pre-packaged human genomes in order to minimize any unexpected errors. However, if your reference genome is not human or  slightly different (i.e. a different hs38 variant) from our pre-packaged ones (and you do not want to re-align), you can make a new reference genome folder by following instructions from https://www.yfish.org/display/PUB/Accucopy#Accucopy-3.7Makeyourownreferencegenomepackage.
