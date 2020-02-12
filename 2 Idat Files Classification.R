library(tidyverse)
dataDirectory <- "~/Desktop/Graduate/MethData/"
# make tumor gland idat dir
tumor_gland_dir <- data.frame(union(paste("mkdir", "tumor_gland"), 
                                    paste("cp", paste0(tumor_gland$sentrix_id, "_",
                                                       tumor_gland$sentrix_position, "*.*"), "tumor_gland")))
write.table(tumor_gland,file="tumor_gland_dir.csv", quote = F, row.names = F, col.names = F)
write.table(tumor_gland_dir,file="tumor_gland_dir.sh", quote = F, row.names = F, col.names = F)
# make colon crypt idat dir
colon_crypt_dir <- data.frame(union(paste("mkdir", "colon_crypt"), 
                                    paste("cp", paste0(colon_crypt$sentrix_id, "_",
                                                       colon_crypt$sentrix_position, "*.*"), "colon_crypt")))
write.table(colon_crypt,file="colon_crypt_dir.csv", quote = F, row.names = F, col.names = F)
write.table(colon_crypt_dir,file="colon_crypt_dir.sh", quote = F, row.names = F, col.names = F)
# make bulk idat dir
bulk_dir <- data.frame(union(paste("mkdir", "bulk"), 
                                    paste("cp", paste0(bulk$sentrix_id, "_",
                                                       bulk$sentrix_position, "*.*"), "bulk")))
write.table(bulk,file="bulk_dir.csv", quote = F, row.names = F, col.names = F)
write.table(bulk_dir,file="bulk_dir.sh", quote = F, row.names = F, col.names = F)
