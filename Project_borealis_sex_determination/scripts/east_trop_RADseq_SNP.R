library(tidyverse)
library(vcfR)

dat <- read_tsv("Project_borealis_sex_determination/data/clustering_analysis/mpile_dp15_tropicalisFamilyEast_sexLinkage.table.gz")

#the filter will give you useful sites, and can then filter on the "sex_chi" column
#the mutate will compute p-values (FDR corrected and not) for association with sex based on the chi values.
useful<-(dat %>% filter( `[1]CHROM` == "Chr07" & 
                  site_status == "goodTag" & 
                  site_type == "hetDad") 
            %>%  mutate(sex_fdr = p.adjust(1-pchisq(sex_chi,1), method = "BH"), 
                        raw_p = 1-pchisq(sex_chi,1)) 
            %>% filter(`[2]POS` == 902118 | 
                       `[2]POS` == 1052539 |  
                       `[2]POS` == 2178835 | 
                       `[2]POS` == 2178837)
        )

library(Biostrings)
a <- readDNAStringSet("Project_borealis_sex_determination/data/genome/XT9_1.fa")
b <- subseq(Chr07, start=901618, stop=902618, width=1000)
writeXStringSet(b, file="Project_borealis_sex_determination/data.fa")