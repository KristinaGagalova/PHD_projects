# Software installation on Genesis

## Install LinuxBrew on genesis

LinuxBrew dependencies [here](http://linuxbrew.sh/)

- Upgrade your gcc version to 4.4 or newer. Follow [this instructions](http://openwall.info/wiki/internal/gcc-local-build). gmp, mpc and mpfr neeed to be installed as well. I did it in my home directory but there are other versions of ggc in the ```/usr/..``` directory. 

- Install Ruby 1.8.6 or newer

- Install cURL. Consider to compile it with *openssl-0.9.8o* or higher because of problems with certificate authentication. In case cURL is not compiled with the correct openssl, you will receive the following error:    
```error:0D0C50A1:asn1 encoding routines:ASN1_item_verify:unknown message digest algorithm``` and you wont be able to download any software packages.       
Check the cURL and openssl versions (**curl --version**).       
    
**WORKAROUND**: Use cURL with *insecure* option. This is not recommended but it works! In order to do so type this command:
```echo insecure >> ~./curlrc```

Once you have installed all the dependencies you are ready to insatall Linuxrew. Now you can update gcc and cURL with Brew.     

## Linuxbrew settings

You may use this option for some formulas but generally avoid to use it.

```
export HOMEBREW_BUILD_FROM_SOURCE=1
```

## Bugs/errors

When installing gcc and its dependencies with brew (ex: glibc) it may complain that there is no gcc installed (chicken or the egg problem). 
In this case add a synlink of gcc in ```~/.linuxbrew/bin``` named gcc-4.4.    
```
ln -s /usr/bin/gcc44 ~/.linuxbrew/bin/gcc-4.4
ln -s /usr/bin/g++44 ~/.linuxbrew/bin/g++-4.4
```

## Install Abyss

The Abyss version installed through LinuxBrew has some incompatibility with genesis openmpi. A recurrent error is:

```
[warn] Epoll ADD(1) on fd 0 failed. Old events were 0; read change was 1 (add); write change was 0 (none); close change was 0 (none): Operation not permitted
[warn] Epoll ADD(4) on fd 1 failed. Old events were 0; read change was 0 (none); write change was 1 (add); close change was 0 (none): Operation not permitted
```

It is recommended to install and compile Abyss independently from brew, following the parameters:

```
#!/bin/bash
#$ -S /bin/bash
#$ -q all.q,mpi.q
#$ -pe ncpus 1
#$ -l mem_token=3.8G,mem_free=3.8G,h_vmem=3.8G
#$ -j y

#code from Austin Hammond

./configure --prefix=$HOME/src/abyss-2.0.2-build/k256 \
    --enable-maxk=256 \ #larger kmer size
    --with-mpi=/usr/mpi/gcc/openmpi-1.4.1/lib64/ \ #use recommended mpi version
    --with-boost=$HOME/src/boost_1_59_0 \ #boost libraries
    CPPFLAGS=-I$HOME/src/sparsehash-build/include/ \ #include sparsehash
    && make && make install
``` 
Remember to set properly the PATHs to this Abyss version.

## Useful links
Brew installation for CentoS [here](https://github.com/Linuxbrew/brew/wiki/CentOS5).     
