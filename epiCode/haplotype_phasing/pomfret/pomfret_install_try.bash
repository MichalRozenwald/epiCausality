# Clone and build Pomfret
cd /home/michalula/software/pomfret_haplotype_phase_ont
screen -S pomfret_haplotype

git clone https://git.oxfordnanolabs.local/xfeng/pomfret/
cd pomfret
make







# 

git clone https://git.oxfordnanolabs.local/xfeng/pomfret/
 git clone https://git.oxfordnanolabs.local/xfeng/pomfret/
Cloning into 'pomfret'...
fatal: unable to access 'https://git.oxfordnanolabs.local/xfeng/pomfret/': Could not resolve host: git.oxfordnanolabs.local

git clone https://github.com/nanoporetech/pomfret.git
cd pomfret

(whatshap-env) michalula@Skynet:~/software/pomfret_haplotype_phase_ont/pomfret$ make
gcc -g -O2 -Wall -Wno-error -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unused-label  -o pomfret kthread.c kstring.c cli.c blockjoin.c main.c -lz -pthread -L /usr/bin/htslib -lhts -Wl,-rpath,'/usr/bin/htslib' -I /usr/bin/htslib
blockjoin.c:7:10: fatal error: htslib/sam.h: No such file or directory  
    7 | #include "htslib/sam.h"
      |          ^~~~~~~~~~~~~~
compilation terminated.

conda activate whatshap-env
sudo apt-get update
sudo apt-get install libhts-dev

=====  Build Instructions
Pomfret builds on Linux and macOS. Dependencies:

gcc
zlib
HTSlib

sudo apt-get update
sudo apt-get install build-essential zlib1g-dev libhts-dev
rm -rf cd /home/michalula/software/pomfret_haplotype_phase_ont/pomfret


# Clone and build Pomfret
cd /home/michalula/software/pomfret_haplotype_phase_ont

git clone https://git.oxfordnanolabs.local/xfeng/pomfret/
git clone https://git.oxfordnanolabs.local/xfeng/pomfret/
Cloning into 'pomfret'...
fatal: unable to access 'https://git.oxfordnanolabs.local/xfeng/pomfret/': Could not resolve host: git.oxfordnanolabs.local

sudo apt-get update
sudo apt-get install build-essential zlib1g-dev libhts-dev

conda activate whatshap-env
conda install -c bioconda htslib zlib

cd /home/michalula/software/pomfret_haplotype_phase_ont
git clone https://github.com/nanoporetech/pomfret.git
cd pomfret
make

   sudo apt-get remove libhts-dev

Clone and build the latest HTSlib:
      cd ~
   git clone https://github.com/samtools/htslib.git
   cd htslib
   git submodule update --init --recursive
   makea
   sudo make install


