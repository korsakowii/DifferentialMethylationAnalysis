# Sex infer
inferSex(ssets_tumor_gland[[1]])
inferSexKaryotypes(ssets_colon_crypt[[1]])
# Ethnicity infer
inferEthnicity(ssets_tumor_gland[[1]])
# Age infer
predictAgeHorvath353(betas_tumor_gland)
# Mean intensity
meanIntensity(ssets_tumor_gland[[1]])
# Bisulfite conversion control using GCT scores
# The closer the score to 1.0, the more complete the bisulfite conversion.
bisConversionControl(ssets_tumor_gland[[1]])
# Visualize all probes from a gene
visualizeGene('DNMT1', betas_tumor_gland, platform='EPIC')
# Visualize probes from arbitrary region
visualizeRegion(
  'chr19',10260000,10380000, betas_tumor_gland, platform='EPIC',
  show.probeNames = FALSE)
# Visualize by probe names
visualizeProbes(c("cg02382400", "cg03738669"), betas_tumor_gland, platform='EPIC')
# Copy number variation
ssets.normal <- sesameDataGet('EPIC.5.normal')
segs <- cnSegmentation(ssets_tumor_gland[[1]], ssets.normal)
segs <- cnSegmentation(ssets_colon_crypt[[1]], ssets.normal)
segs <- cnSegmentation(ssets_tumor_gland[[1]], ssets_colon_crypt)
visualizeSegments(segs)

