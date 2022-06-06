#!/usr/bin/bash

#Based on: https://stackoverflow.com/questions/41842304/reshaping-from-wide-to-long-format

'
Input format, 2 columns
A1  A4|A5|A6  
B1  B4|B5|B6  
C1  C4|C5|C6

Output format
A1  A4
A1  A5
A1  A6  
B1  B4
B1  B5
B1  B6  
C1  C4
C1  C5
C1  C6

Deliniter is "|"
'

input_file=$1

while read -r f1 c1; do
  # split the comma delimited field 'c1' into its constituents
  for c in ${c1//|/ }; do
     printf "$f1\t$c\n"
  done
done < $input_file
