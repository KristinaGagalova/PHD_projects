#!/usr/bin/bash

go_2columns=$1
name=$2
 
awk -v NAM=$name '{print NAM,$1,$1,"optional",$2,"pubmed","unknown",$2,"unknown","optional","optional","protein","unknown",NAM,$1,"optional","optional"}' $go_2columns | tr ' ' '\t'
