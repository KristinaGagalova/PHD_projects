## Install Perl modules on filer

Through CPAN

[Gin Link](http://gin.bcgsc.ca/plone/groups/lims/Wiki/Standards%20docs/Installing%20Perl%20Modules/)

```
perl -MCPAN -e 'install HTML::Template'
```

```
#example
cpan install XML::Simple
```

cpan is configured to provide one setting to Makefile.PL, this may create a conflict with local::lib, which uses the INSTALL_BASE paradigm. Use one or the other      

How to stop cpan of using INSTALL_BASE:

```
o conf makepl_arg ''
o conf mbuildpl_arg ''
o conf commit
```

## Eventual errors

### Error in Perl versions

```
perl: symbol lookup error: /home/kgagalova/perl5/lib/perl5/x86_64-linux-thread-multi/auto/Cwd/Cwd.so: undefined symbol: Perl_Istack_sp_ptr
```

Major Perl versions (5.6, 5.8, 5.10, 5.12) aren't binary compatible with each other. You can't use a module with a compiled component if you're using a different major version of Perl than the one that was used to install the module. (I think you can also have problems if you use an earlier minor version of Perl than the one with which the module was installed.)

### See where perl libraries are installed

```
perldoc -l Time::HiRes
```

### install independent version of Perl

```
wget http://www.cpan.org/src/5.0/perl-5.20.1.tar.gz
tar -xzf perl-5.20.1.tar.gz
cd perl-5.20.1
./Configure -des -Dprefix=$HOME/localperl -Dusethreads #add multithreading
make
make test
make install
export PATH=$HOME/localperl/bin:$PATH
```

```
curl -L http://cpanmin.us | perl - App::cpanminus
cpanm Module::Name
```


## install modules in specif directory with CPAN

```
#to start into the local directory installation
PERL5LIB=/home/kgagalova/localperl/lib/site_perl/5.20.1/x86_64-linux-thread-multi perl -MCPAN -e shell
```

### Perl modules

Check the link [here](https://stackoverflow.com/questions/2526804/how-is-perls-inc-constructed-aka-what-are-all-the-ways-of-affecting-where-pe)

## Current setting on filer

Perl modules
```
PERL5LIB=/home/kgagalova/localperl/lib/site_perl/5.20.1/x86_64-linux-thread-multi
```

Perl version

```
PATH=$HOME/localperl/bin:$PATH && perl run  
```

Config args, include threading

```
config_args='-des -Dprefix=/home/kgagalova/localperl -Dusethreads -Dotherlibdirs=/home/kgagalova/localperl/lib/site_perl/5.20.1/x86_64-linux-thread-multi'
```
