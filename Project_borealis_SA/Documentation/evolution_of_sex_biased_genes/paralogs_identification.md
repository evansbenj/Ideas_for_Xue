# Identify paralogs 
For the phylogenetic analysis, I would need to identify paralogs in *X. laevis* and and in *X. borealis*. 

*Collasped transcripts
Collapsed transcripts will be mapped againsted its own using blastn. The reason of not using Gmap is because Gmap only report the top alignment (and chimeric) but for paralogs, we the the second chr/scaffold that the transcripts align to and not the first one. 
```
#run the blastn, pipe the blastn output to the perl script
```
