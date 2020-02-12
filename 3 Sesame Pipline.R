library(sesame)
data_dir <- "~/Desktop/Graduate/MethData/"
# Read IDATs into SigSet list
ssets_tumor_gland <- lapply(
  searchIDATprefixes(paste0(data_dir, "tumor_gland")),
  readIDATpair)
ssets_colon_crypt <- lapply(
  searchIDATprefixes(paste0(data_dir, "colon_crypt")),
  readIDATpair)
# Background subtraction
# normal-exponential deconvolution using out-of-band probes `noob`
# `noobsb` regressing out the green-to-red and red-to-green relationship using Type-I probes
sset_tumor_gland.nb <- ssets_tumor_gland[[1]] %>% noob() %>% noobsb() 
sset_colon_crypt.nb <- ssets_colon_crypt[[1]] %>% noob() %>% noobsb()

# Type-I channel inference
sset_tumor_gland.nb.TypeICorrected <- inferTypeIChannel(sset_tumor_gland.nb)
sset_colon_crypt.nb.TypeICorrected <- inferTypeIChannel(sset_colon_crypt.nb)

# Dye bias correction
# Residual dye bias can be corrected using nonlinear quantile interpolation with Type-I probes.
sset_tumor_gland.dbNlinear <- dyeBiasCorrTypeINorm(sset_tumor_gland.nb.TypeICorrected)
qqplot(
  slot(sset_tumor_gland.dbNlinear, 'IR'), slot(sset_tumor_gland.dbNlinear, 'IG'),
  xlab='Type-I Red Signal', ylab='Type-I Grn Signal',
  main='Nonlinear Correction', cex=0.5)
abline(0,1,lty='dashed')

sset_colon_crypt.dbNlinear <- dyeBiasCorrTypeINorm(sset_colon_crypt.nb.TypeICorrected)
qqplot(
  slot(sset_colon_crypt.dbNlinear, 'IR'), slot(sset_colon_crypt.dbNlinear, 'IG'),
  xlab='Type-I Red Signal', ylab='Type-I Grn Signal',
  main='Nonlinear Correction', cex=0.5)
abline(0,1,lty='dashed')

# Get betas
betas_tumor_gland <- getBetas(sset_tumor_gland.dbNlinear)
betas_colon_crypt <- getBetas(sset_colon_crypt.dbNlinear)
# Extra SNP allele frequencies
extraSNPAFs_tumor_gland <- getAFTypeIbySumAlleles(sset_tumor_gland.dbNlinear)
extraSNPAFs_colon_crypt <- getAFTypeIbySumAlleles(sset_colon_crypt.dbNlinear)