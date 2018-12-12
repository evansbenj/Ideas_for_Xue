# Collapsing transcripts using corset
Corset clusters transcripts from de novo transcriptome assembly to gene-level counts based on shared reads and expression patterns. 
- Corset documentation: https://github.com/Oshlack/Corset/wiki/InstallingRunningUsage
- Example: https://github.com/Oshlack/Corset/wiki/Example

### salmon
building index
```bash
salmon index --index Trinity --transcripts trinity_out_dir.Trinity.fasta
```
quantifying reads
```bash
##########bash script to run read quantificaitons for each individual separately using salmon  
#!/bin/bash

FILES=`ls /home/benf/Borealis-Family-Transcriptomes-July2017/Data/Trimmed/BJE*_liver_R1_scythe.fastq.gz | sed 's/_liver_R1_scythe.fastq.gz//g'`
for F in $FILES ; do
        R1=${F}_liver_R1_scythe.fastq.gz
        R2=${F}_Rliver_R1_scythe.fastq.gz
        salmon quant --index Trinity --libType IU --dumpEq -1 $R1 -2 $R2 --output ${F}.out
done
```
### corset
```bash
corset -g 1,1,1,2,2,2 -n A1,A2,A3,B1,B2,B3 -i salmon_eq_classes BJE*/aux_info/eq_classes.txt
```
```bash
for FILE in `ls *.bam` ; do
   corset -r true-stop $FILE &
done
wait
```
