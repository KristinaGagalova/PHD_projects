#Description: http://avrilomics.blogspot.com/2015/06/creating-go-slim-and-mapping-go-terms.html

input_gaf=$1

owlTools=/home/kgagalova/src/owltools

oboBasic=/projects/spruceup_scratch/dev/SprucePaper2018/GoTermsAnalysis/MapToGOslim/Dependencies/go-basic.obo
goSlim=/projects/spruceup_scratch/dev/SprucePaper2018/GoTermsAnalysis/MapToGOslim/Dependencies/goslim_plant_terms.txt

#test
java -jar $owlTools $oboBasic --gaf $input_gaf --map2slim \
	--idfile $goSlim --write-gaf my_Slimgo_terms_plant.gaf

