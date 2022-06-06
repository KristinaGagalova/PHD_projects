## Run Maker for gene annotation of Kollector output

Since a large part of the Maker genes are partial and probably some starting/ending exons are missing, it is worth to try to recover some of the genes by running Maker on the reconstructed by Kollector contigs.   


1) Run RepeatMasker - Use personalized libraries from Austin's LTR-retriever and other tools.     
Locations libraries:    
```
##WS
/projects/spruceup/scratch/pglauca/WS77111/annotation/repeat-elements/WS77111-combined-repeatLib-plusRepBase22.08.fa
##PG29 - same used for Q903
/projects/spruceup/scratch/interior_spruce/PG29/annotation/repeat-elements/PG29-combined-repeatLib-plusRepBase22.08.fa
```

It is recommended to chop the file into batches since RepeatMasker is quite slow.      
Also: fasta headers can be max 50 characters. You may trim them or rename (keep always file with original names!!)    

Independent run of RepeatMasker before Maker
```
#lib spruce_rep is the customized library, fasta is the batch of fasta file
RepeatMasker -pa 12 -e ncbi -lib $spruce_rep -dir $Nam $fasta
```

Reformat output

```
# create GFF3
rmOutToGFF3.pl Full_mask/Boa_constrictor_SGA_7C.scaffolds.full_mask.out > Full_mask/Boa_constrictor_SGA_7C.scaffolds.full_mask.out.gff3
# isolate complex repeats
grep -v -e "Satellite" -e ")n" -e "-rich" Boa_constrictor_SGA_7C_scaffolds.full_mask.gff3 \
  > Boa_constrictor_SGA_7C_scaffolds.full_mask.complex.gff3
# reformat to work with MAKER
cat Boa_constrictor_SGA_7C_scaffolds.full_mask.complex.gff3 | \
  perl -ane '$id; if(!/^\#/){@F = split(/\t/, $_); chomp $F[-1];$id++; $F[-1] .= "\;ID=$id"; $_ = join("\t", @F)."\n"} print $_' \
  > Boa_constrictor_SGA_7C_scaffolds.full_mask.complex.reformat.gff3
```

2) Run Maker

Create correctly formatted Maker files     
```
maker -CTL
```

This will output 

- maker_bopts.ctl - Blast options
- maker_opts.ctl - Maker input and options
- maker_exe.ctl - maker executables


