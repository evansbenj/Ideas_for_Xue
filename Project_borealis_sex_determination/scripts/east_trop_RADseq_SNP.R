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

#library(Biostrings)
#a <- readDNAStringSet("Project_borealis_sex_determination/data/genome/XT9_1.fa")
#b <- subseq(Chr07, start=901618, stop=902618, width=1000)
#writeXStringSet(b, file="Project_borealis_sex_determination/data.fa")

male = c("[23]BJE4689_boy_GhE.fq.gz:GT", "[25]BJE4691_boy_GhE.fq.gz:GT", 
         "[27]BJE4693_boy_GhE.fq.gz:GT","[29]BJE4697_boy_GhE.fq.gz:GT")
female = c("[18]BJE4684_girl_GhE.fq.gz:GT","[19]BJE4685_girl_GhE.fq.gz:GT", 
           "[20]BJE4686_girl_GhE.fq.gz:GT","[21]BJE4687_girl_GhE.fq.gz:GT", 
           "[22]BJE4688_boy_GhE.fq.gz:GT", "[24]BJE4690_girl_GhE.fq.gz:GT",
           "[26]BJE4692_girl_GhE.fq.gz:GT", "[28]BJE4694_girl_GhE.fq.gz:GT")

#[31]Xt_mom_BJE4361:GT
#[30]Xt_dad_BJE4362:GT

additional <-(dat %>% filter( site_type == "hetDad"
                              , `[1]CHROM` == "Chr07" 
                              , site_status == "goodTag" 
                              ) 
              %>% mutate(
                          male_match = (`[30]Xt_dad_BJE4362:GT` == `[23]BJE4689_boy_GhE.fq.gz:GT`) 
                                      + (`[30]Xt_dad_BJE4362:GT` == `[25]BJE4691_boy_GhE.fq.gz:GT`)
                                      + (`[30]Xt_dad_BJE4362:GT` == `[27]BJE4693_boy_GhE.fq.gz:GT`)
                                      + (`[30]Xt_dad_BJE4362:GT` == `[29]BJE4697_boy_GhE.fq.gz:GT`)
                         , female_match = (`[31]Xt_mom_BJE4361:GT` == `[18]BJE4684_girl_GhE.fq.gz:GT` | `[18]BJE4684_girl_GhE.fq.gz:GT` == "./.")
                                        + (`[31]Xt_mom_BJE4361:GT` == `[19]BJE4685_girl_GhE.fq.gz:GT` | `[19]BJE4685_girl_GhE.fq.gz:GT` == "./." )
                                        + (`[31]Xt_mom_BJE4361:GT` == `[20]BJE4686_girl_GhE.fq.gz:GT` | `[20]BJE4686_girl_GhE.fq.gz:GT` == "./.")
                                        + (`[31]Xt_mom_BJE4361:GT` == `[21]BJE4687_girl_GhE.fq.gz:GT` | `[21]BJE4687_girl_GhE.fq.gz:GT`== "./.")
                                        + (`[31]Xt_mom_BJE4361:GT` == `[22]BJE4688_boy_GhE.fq.gz:GT` | `[22]BJE4688_boy_GhE.fq.gz:GT` == "./.")
                                        + (`[31]Xt_mom_BJE4361:GT` == `[24]BJE4690_girl_GhE.fq.gz:GT`| `[24]BJE4690_girl_GhE.fq.gz:GT` == "./.")
                                        + (`[31]Xt_mom_BJE4361:GT` == `[26]BJE4692_girl_GhE.fq.gz:GT`| `[26]BJE4692_girl_GhE.fq.gz:GT` == "./.")
                                        + (`[31]Xt_mom_BJE4361:GT` == `[28]BJE4694_girl_GhE.fq.gz:GT`| `[28]BJE4694_girl_GhE.fq.gz:GT` == "./.")
                         , father_daughter_match = (`[30]Xt_dad_BJE4362:GT` == `[18]BJE4684_girl_GhE.fq.gz:GT` )
                                         + (`[30]Xt_dad_BJE4362:GT` == `[19]BJE4685_girl_GhE.fq.gz:GT`  )
                                         + (`[30]Xt_dad_BJE4362:GT` == `[20]BJE4686_girl_GhE.fq.gz:GT` )
                                         + (`[30]Xt_dad_BJE4362:GT` == `[21]BJE4687_girl_GhE.fq.gz:GT` )
                                         + (`[30]Xt_dad_BJE4362:GT` == `[22]BJE4688_boy_GhE.fq.gz:GT` )
                                         + (`[30]Xt_dad_BJE4362:GT` == `[24]BJE4690_girl_GhE.fq.gz:GT`)
                                         + (`[30]Xt_dad_BJE4362:GT` == `[26]BJE4692_girl_GhE.fq.gz:GT`)
                                         + (`[30]Xt_dad_BJE4362:GT` == `[28]BJE4694_girl_GhE.fq.gz:GT`)
                         , total_match = male_match +female_match
                         , sex_fdr = p.adjust(1-pchisq(sex_chi,1), method = "BH")
                         , raw_p = 1-pchisq(sex_chi,1)
                        )
              %>% filter( #male_match >=1
                          #,female_match >= 7
                          #,total_match >= 8
                          father_daughter_match <=0
                          ,raw_p <=0.5
                        )
              %>% select (`[1]CHROM`, `[2]POS`,female_match, male_match, total_match,father_daughter_match, raw_p,`[3]REF` 
                          , male_4362_dad =`[30]Xt_dad_BJE4362:GT`
                          , male_4691 = `[25]BJE4691_boy_GhE.fq.gz:GT`
                          , male_4693 = `[27]BJE4693_boy_GhE.fq.gz:GT`
                          , male_4697 = `[29]BJE4697_boy_GhE.fq.gz:GT`
                          , male_4689 = `[23]BJE4689_boy_GhE.fq.gz:GT`
                          , female_4361_mom = `[31]Xt_mom_BJE4361:GT`
                          , female_4684 = `[18]BJE4684_girl_GhE.fq.gz:GT`
                          , female_4685 = `[19]BJE4685_girl_GhE.fq.gz:GT`
                          , female_4686 = `[20]BJE4686_girl_GhE.fq.gz:GT`
                          , female_4687 = `[21]BJE4687_girl_GhE.fq.gz:GT`
                          , female_4688 = `[22]BJE4688_boy_GhE.fq.gz:GT`
                          , female_4690 = `[24]BJE4690_girl_GhE.fq.gz:GT`
                          , female_4692 = `[26]BJE4692_girl_GhE.fq.gz:GT`
                          , female_4694 = `[28]BJE4694_girl_GhE.fq.gz:GT`
                          )
              
)
