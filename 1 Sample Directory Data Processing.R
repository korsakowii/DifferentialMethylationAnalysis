library(tidyverse)
library(unpivotr)
dataDirectory <- "~/Desktop/Graduate/MethData/"
directory <- read_csv(paste0(dataDirectory,"directory.csv"))
# Actually we can delete the column "type"
# SI == small intestine
# duo == Duodenum
# crypt(stem cells)
# G2P2 == two pregnancies and two deliveries
# tumorSURF == tumor surface?
directory <- directory %>% select(-type)
unique(directory$description)

# tumor gland sample
tumor_gland <- directory %>% filter(description == "tumor gland")
## count # of sample / person
table(str_sub(tumor_gland$sample_id, end = 1))
## count # of person for sample
table(str_sub(tumor_gland$sample_id,start = 2))


# normal tissue
colon_crypt <- directory %>%
  filter((grepl("colon", tissue)&grepl("crypt", description)))
## count # of sample / person
table(str_sub(colon_crypt$sample_id, end = 1))
## count # of person for sample
table(str_sub(colon_crypt$sample_id,start = 2))

small_crypt <- directory %>% 
  filter(grepl("(SI crypt|crypt SI)", description))
endom_gland <- directory %>% 
  filter(grepl("(endo|endometrial)", description))

# bulk sample
bulk_sample <- directory %>% filter(description == "bulk")
## count # of sample / person
table(str_sub(bulk_sample$sample_id, end = 1))
## count # of person for sample
table(str_sub(bulk_sample$sample_id,start = 2))

bulk_group <- bulk_sample  %>% arrange(reverse(sample_id))

write.table(bulk_group, file = "bulk_group.csv", quote = F, sep = ",", row.names = F)


