## Install BBT

```
git clone git@github.com:bcgsc/biobloom.git
```
Use -b for specific branches

```
git submodule update --init
```

Use libraries from Justin here or linuxbrew
```
export PATH=$PATH:/gsc/btl/linuxbrew/bin:/home/cjustin/arch/x86-64/bin/
```

```
.autogen.sh
./configure --prefix=/path/here
make
make install
```


