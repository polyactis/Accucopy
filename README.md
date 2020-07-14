# Accucopy

Accucopy is a CNA-calling method that extends our previous Accurity model to predict both total (TCN) and allele-specific copy numbers (ASCN) for the tumor genome. Accucopy adopts a tiered Gaussian mixture model coupled with an innovative autocorrelation-guided EM algorithm to find the optimal solution quickly. The Accucopy model utilizes information from both total sequencing coverage and allelic sequencing coverage. Through comparative analyses in both simulation and real-sequencing samples, we demonstrate that Accucopy is more accurate than existing methods

The manuscript, titled "Accucopy: accurate and fast inference of allele-specific copy number alterations from low-coverage low-purity tumor sequencing data", is under review.


# License

The license follows our institute policy that you can use the program for free as long as you are using Accucopy strictly for non-profit research purposes. However, if you plan to use Accucopy for commercial purposes, a license is required and please contact yuhuang@simm.ac.cn or polyactis@gmail.com to obtain one.

The full-text of the license is included in the software package.

# Compile

To compile, go into src_o/ and type "make debug" or "make release".

Accucopy is not easy to compile. If you just want to use Accucopy, please use its docker, https://www.yfish.org/display/PUB/Accucopy#Accucopy-Docker.

For more information, check https://www.yfish.org/display/PUB/Accucopy.
