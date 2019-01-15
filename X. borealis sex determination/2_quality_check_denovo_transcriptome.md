# *de novo* Assembled Transcriptome Quality Check
## Transfer Data to Info
```bash
#I transfer the transcriptome and supertranscriptome data, and mapping result back to Info for down stream analysis
#since on Info, jobs will be ran more immediately without the queueing system, which Graham use and sometime jobs just queueing for a day or more.
scp -r songxy@graham.computecanada.ca:/home/songxy/projects/def-ben/songxy/tropicalis_gonad_transcriptome/data/tropicalis_gonad_supertranscriptome_dec2018/ /home/xue/tropicalis_gonad_transcriptome_Dec2018/data/tropicali_gonad_transcriptome_trinityOut

scp -r songxy@graham.computecanada.ca:/home/songxy/projects/def-ben/songxy/tropicalis_gonad_transcriptome/data/tropicalis_transcriptome_build_dec2018 /home/xue/tropicalis_gonad_transcriptome_Dec2018/data/tropicali_gonad_transcriptome_trinityOut
```

## basic check
I ran the below script provided by Trinity to do the quality check
```
```
