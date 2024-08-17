setwd("C:/Users/Atefeh/Desktop/ovarian")
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)



lnames1 = load(file = "Stroma.RData");
lnames1


lnames2 = load(file = "Tumor.RData");
lnames2




lnamesStroma=load(file = "Stromamodules.RData")
lnamesStroma


lnamesTumor=load(file = "Tumormodules.RData")
lnamesTumor

setLabels = c("StromaN", "TomurN");
multiExpr = list( StromaN = list(data = data.Expression.Stroma),TumorN = list(data = data.Expression.Tumor));
multiColor = list(StromaN = moduleColors1);




system.time( {
  mp = modulePreservation(multiExpr, multiColor,
                          referenceNetworks = 1,
                          nPermutations =50,
                          randomSeed = 1,
                          quickCor = 0,
                          verbose = 3)
} );
# Save the results
save(mp, file = "modulePreservation.RData");

