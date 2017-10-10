Purpose: design primer to amplify the maybe sex-linked Borealis transcript

# Primer Blast
I used the scafold sequence as input query for Primer Blast. Specified that the forward primer range as 961-1049 and reverse primer range as 1400-1509.
The top 10 primers were shown. I went through them manually to see which one was the best based on the following criteria:
  - didnt't start with long "TTTT"
  - have lower number of self-completarity
  - melting point closer to 60C
  - reasonable GC content

The below primer pair is the best among the result after consider the above criteria:
```
TAGTTACAAATGGTGGGTTCTGA
CACTTTGCCTCACCCTCATAG

```
with the following statistic:
```
                 Sequence (5'->3')	Template strand_Length	Start	Stop	Tm	GC%	Self_complementarity	Self_3'_complementarity
Forward primer	TAGTTACAAATGGTGGGTTCTGA	Plus	23	997	1019	57.43	39.13	3.00	3.00
Reverse primer	CACTTTGCCTCACCCTCATAG	Minus	21	1456	1436	58.36	52.38	3.00	2.00
Product length	460
```
The above was chosen among the others(below):
```

Sequence (5'->3')	Template strand	Length	Start	Stop	Tm	GC%	Self complementarity	Self 3' complementarity
Forward primer	ACTAGTTACAAATGGTGGGTTCT	Plus	23	995	1017	57.42	39.13	6.00	0.00
Reverse primer	CACTTTGCCTCACCCTCATA	Minus	20	1456	1437	57.21	50.00	2.00	2.00
Product length	462

Sequence (5'->3')	Template strand	Length	Start	Stop	Tm	GC%	Self complementarity	Self 3' complementarity
Forward primer	GCATGAGTATGCATCTATTCCTGT	Plus	24	970	993	58.76	41.67	6.00	2.00
Reverse primer	AGATTCAGGCCCAGCACTTT	Minus	20	1470	1451	59.59	50.00	4.00	2.00
Product length	501

Sequence (5'->3')	Template strand	Length	Start	Stop	Tm	GC%	Self complementarity	Self 3' complementarity
Forward primer	TCACTAGTTACAAATGGTGGGTT	Plus	23	993	1015	57.69	39.13	6.00	0.00
Reverse primer	GCACTTTGCCTCACCCTCATA	Minus	21	1457	1437	60.34	52.38	3.00	2.00
Product length	465

Sequence (5'->3')	Template strand	Length	Start	Stop	Tm	GC%	Self complementarity	Self 3' complementarity
Forward primer	TGAGTATGCATCTATTCCTGTCAC	Plus	24	973	996	58.15	41.67	6.00	3.00
Reverse primer	GAGATTCAGGCCCAGCACTTT	Minus	21	1471	1451	60.62	52.38	4.00	1.00
Product length	499
```
