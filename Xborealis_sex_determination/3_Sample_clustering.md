# Samples clustering - separate male and female samples
## Principle Component Analysis (PCA)
### Learning
I know nothing about PCA so I started this analysis by learning what it is and how to do it in R. 

#### What is PCA?

Reference
- StatQuest explained it very well, so I recommended watching this first:https://www.youtube.com/watch?v=FgakZw6K1QQ  
- Bio-Star Q&A: https://www.biostars.org/p/171925/
- nature methods: https://www.nature.com/articles/nmeth.4346
- Stackexchange Q&A: https://stats.stackexchange.com/questions/2691/making-sense-of-principal-component-analysis-eigenvectors-eigenvalues
- https://www.biostars.org/p/280615/


#### How to do PCA?
In the paper by Fehrmann et al (2015), their "preprocessing and aggregation of raw data were performed according to the robust multi-array average (RMA) algorithm in combination with quantile normalization".

Reference
- StatQuest explain the basic very well, I recommended watching this first: https://statquest.org/2017/11/27/statquest-pca-in-r-clearly-explained/
- https://www.nature.com/articles/ng.3173#methods
- https://www.hindawi.com/journals/bmri/2018/2906292/
- https://www.datacamp.com/community/tutorials/pca-analysis-r
- https://www.biostars.org/p/282685/
- http://strata.uga.edu/8370/handouts/pcaTutorial.pdf
- http://uc-r.github.io/pca


## Other method
PCA analysis didn't really work in seperating our samples into two clusters. We are going to try more unsupervised clustering methods. A biostars Q&A (https://www.biostars.org/p/56897/) provided good resource to check out. 
### kmeans
A phd student from Dr. Brian Golding's lab suggested kmean and set k=2 if we are looking for two group of clustering (males and females). 

### By examining dmrt1 expression

```
#blastn index tropicalis de novo transcriptome 


# blastn dmrt1 genes into tropicalis transcriptome


```
