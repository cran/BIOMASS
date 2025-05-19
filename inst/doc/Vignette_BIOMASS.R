## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE, echo = TRUE,
  comment = "#>", fig.align = "center")
require(knitr)
require(BIOMASS)

## ----load_NouraguesTrees------------------------------------------------------
data("NouraguesTrees")
knitr::kable(head(NouraguesTrees))

## ----eval=FALSE---------------------------------------------------------------
#  # By default
#  createCache()
#  # Or if you want to set your own cache folder
#  createCache("the_path_to_your_cache_folder")
#  # Or
#  options("BIOMASS.cache" = "the_path_to_your_cache_folder")
#  

## ----save_correctTaxo, eval=FALSE---------------------------------------------
#  Taxo <- correctTaxo(
#    genus = NouraguesTrees$Genus, # genus also accepts the whole species name (genus + species) or (genus + species + author)
#    species = NouraguesTrees$Species,
#    useCache = TRUE, verbose = FALSE)
#  saveRDS(Taxo, file = "saved_data/Taxo_vignette.rds")

## ----load_correctTaxo, include=FALSE------------------------------------------
Taxo <- readRDS(file = "saved_data/Taxo_vignette.rds")

## ----corrected_taxo-----------------------------------------------------------
NouraguesTrees$GenusCorrected <- Taxo$genusCorrected
NouraguesTrees$SpeciesCorrected <- Taxo$speciesCorrected

## ----species_correction_example-----------------------------------------------
NouraguesTrees$Species[4]
Taxo[4,]

## ----getTaxonomy, eval=FALSE--------------------------------------------------
#  APG <- getTaxonomy(NouraguesTrees$GenusCorrected, findOrder = TRUE)
#  NouraguesTrees$familyAPG <- APG$family
#  NouraguesTrees$orderAPG <- APG$order

## ----getWoodDensity-----------------------------------------------------------
wood_densities <- getWoodDensity(
  genus = NouraguesTrees$GenusCorrected,
  species = NouraguesTrees$SpeciesCorrected,
  stand = NouraguesTrees$Plot # for unidentified or non-documented trees in the reference database
)

NouraguesTrees$WD <- wood_densities$meanWD

## ----wd_informations----------------------------------------------------------
# At species level
sum(wood_densities$levelWD == "species")
# At genus level
sum(wood_densities$levelWD == "genus")
# At plot level
sum(!wood_densities$levelWD %in% c("genus", "species"))

## ----LocalWoodDensity, eval=FALSE---------------------------------------------
#  LocalWoodDensity <- data.frame(
#    genus = c("Paloue", "Handroanthus"),
#    species = c("princeps", "serratifolius"),
#    wd = c(0.65, 0.72) )
#  
#  add_wood_densities <- getWoodDensity(
#    genus = NouraguesTrees$GenusCorrected,
#    species = NouraguesTrees$SpeciesCorrected,
#    family = NouraguesTrees$familyAPG,
#    stand = NouraguesTrees$Plot,
#    addWoodDensityData = LocalWoodDensity
#  )

## ----load_NouraguesHD---------------------------------------------------------
data("NouraguesHD")

## ----multiple_modelHD, message=FALSE------------------------------------------
HD_res <- modelHD(
  D = NouraguesHD$D, H = NouraguesHD$H,
  useWeight = TRUE, drawGraph = T)

kable(HD_res)

## ----log2_modelHD-------------------------------------------------------------
HDmodel <- modelHD(
  D = NouraguesHD$D, H = NouraguesHD$H,
  method = "log2", useWeight = TRUE)

H_model <- retrieveH(
  D = NouraguesTrees$D,
  model = HDmodel)

NouraguesTrees$H_model <-  H_model$H

## ----retrieveH_feldspausch----------------------------------------------------
H_feldspausch <- retrieveH(
  D = NouraguesTrees$D,
  region = "GuianaShield")

NouraguesTrees$H_feldspausch <- H_feldspausch$H

## ----retrieveH_chave, eval=FALSE----------------------------------------------
#  data("NouraguesCoords") #contains corner coordinates
#  coords <- apply(NouraguesCoords[c("Long","Lat")] , 2, mean) # compute the mean of the corner coordinates
#  
#  H_chave  <- retrieveH(
#    D = NouraguesTrees$D,
#    coord = coords)
#  
#  NouraguesTrees$H_chave <- H_chave$H

## ----computeAGB---------------------------------------------------------------
NouraguesTrees$AGB <- computeAGB(
  D = NouraguesTrees$D,
  WD = NouraguesTrees$WD,
  H = NouraguesTrees$H_model #here with the local H-D predictions
  )

## ----save_AGB, include=FALSE--------------------------------------------------
#saveRDS(NouraguesTrees, file = "saved_data/NouraguesTreesAGB.rds")

## ----computAGB_chave, eval=FALSE----------------------------------------------
#  NouraguesTrees$AGB_Chave <- computeAGB(
#      D = NouraguesTrees$D,
#      WD = NouraguesTrees$WD,
#      coord = coords)

## ----D_error_prop-------------------------------------------------------------
D_error_prop <- AGBmonteCarlo(
  D = NouraguesTrees$D, WD = NouraguesTrees$WD, H = NouraguesTrees$H_model,
  Dpropag = "chave2004",
  errWD = rep(0,nrow(NouraguesTrees)), errH = 0 # no error propagation on WD and H here
)

## ----WD_error_prop------------------------------------------------------------
WD_error_prop <- AGBmonteCarlo(
  D = NouraguesTrees$D, WD = NouraguesTrees$WD, H = NouraguesTrees$H_model,
  errWD = wood_densities$sdWD,
  Dpropag = 0 , errH = 0 # no error propagation on D and H here
)

## ----H_model_error_prop-------------------------------------------------------
H_model_error_prop <- AGBmonteCarlo(
  D = NouraguesTrees$D, WD = NouraguesTrees$WD, # we do not provide H
  HDmodel = HDmodel, # but we provide HDmodel
  Dpropag = 0 , errWD = rep(0,nrow(NouraguesTrees)) # no error propagation on D and WD here
)

## ----H_feld_error_prop, eval = FALSE------------------------------------------
#  H_feld_error_prop <- AGBmonteCarlo(
#    D = NouraguesTrees$D, WD = NouraguesTrees$WD,
#    H = NouraguesTrees$H_feldspausch, errH = H_feldspausch$RSE, # we provide H and errH
#    Dpropag = 0 , errWD = rep(0,nrow(NouraguesTrees)) # no error propagation on D and WD here
#  )

## ----H_chave_error_prop, eval = FALSE-----------------------------------------
#  H_chave_error_prop <- AGBmonteCarlo(
#    D = NouraguesTrees$D, WD = NouraguesTrees$WD, # we do not provide H
#    coord = coords, # but we provide the vector of median coordinates of the plots
#    Dpropag = 0 , errWD = rep(0,nrow(NouraguesTrees)) # no error propagation on D and WD here
#  )

## ----error_prop---------------------------------------------------------------
error_prop <- AGBmonteCarlo(
  D = NouraguesTrees$D, WD = NouraguesTrees$WD, # we do not provide H
  HDmodel = HDmodel, # but we provide HDmodel
  Dpropag = "chave2004",
  errWD = wood_densities$sdWD)

error_prop[(1:4)]

## ----summaryByPlot------------------------------------------------------------
AGB_by_plot <- summaryByPlot(AGB_val = error_prop$AGB_simu, plot = NouraguesTrees$Plot, drawPlot = TRUE)


## ----tricks_Hmix, eval = FALSE------------------------------------------------
#  # NouraguesHD contains 163 trees that were not measured
#  NouraguesHD$Hmix <- NouraguesHD$H
#  NouraguesHD$RSEmix <- 0.5
#  filt <- is.na(NouraguesHD$Hmix)
#  NouraguesHD$Hmix[filt] <- retrieveH(NouraguesHD$D, model = HDmodel)$H[filt]
#  NouraguesHD$RSEmix[filt] <- HDmodel$RSE

## ----tricks_Hmix_prop, eval=FALSE---------------------------------------------
#  wd <- getWoodDensity(NouraguesHD$genus, NouraguesHD$species)
#  resultMC <- AGBmonteCarlo(
#    D = NouraguesHD$D, WD = wd$meanWD, errWD = wd$sdWD,
#    H = NouraguesHD$Hmix, errH = NouraguesHD$RSEmix,
#    Dpropag = "chave2004"
#  )
#  summaryByPlot(AGB_val = resultMC$AGB_simu, plot = NouraguesHD$plotId, drawPlot = TRUE)

