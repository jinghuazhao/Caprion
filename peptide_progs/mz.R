#!/usr/local/Cluster-Apps/ceuadmin/R/4.4.1-icelake/bin/Rscript

options(width=200)

# ZWK .raw data
library(rawrr)
spectra_ZWK <- "~/Caprion/pre_qc_data/spectral_library_ZWK"
raw_files <- list.files(spectra_ZWK, pattern = "\\.raw$", full.names = TRUE)
## collectively
suppressMessages(library(MsBackendRawFileReader))
ZWK <- Spectra::backendInitialize(MsBackendRawFileReader::MsBackendRawFileReader(),
       files = raw_files)
class(ZWK)
methods(class=class(ZWK))
Spectra(ZWK)
spectraData(ZWK)
ZWK
ZWK |> Spectra::spectraVariables()
save(ZWK,file="~/Caprion/analysis/work/ZWK.rda")

# ZYQ/UDP
library(utils)
spectra <- "~/Caprion/pre_qc_data/spectra"
zip_files <- dir(spectra, recursive = TRUE, full.names=TRUE)
work_dir <- "~/Caprion/analysis/work"
for (zip_file in zip_files) unzip(zip_file, exdir=work_dir)
ZYQ_UDP <- Spectra::backendInitialize(MsBackendRawFileReader::MsBackendRawFileReader(),
           files = dir(work_dir,patt="raw",full.names=TRUE))
class(ZYQ_UDP)
ZYQ_UDP
ZYQ_UDP |> Spectra::spectraVariables()
save(ZYQ_UDP,file="~/Caprion/analysis/work/ZYQ_UDP.rda")

# mzML
library(mzR)
d <- "/rds/project/rds-MkfvQMuSUxk/interval/caprion_proteomics/spectra"
mz <-  openMSfile(file.path(d,"budp082710c21101xms1.mzML"))
methods(class="mzRpwiz")
spec <- spectra(mz)
class(spec)
length(spec)
lapply(spec,head,3)

# MSnbase
library(MSnbase)
mzXML <- readMSData("szwk901104i19801xms1.mzXML")
extractSpectraData(mzXML)
hasSpectra("szwk901104i19801xms1.mzML,gz")
hasChromatograms("szwk901104i19801xms1.mzML,gz")
plot2d(mzXML,z="peaks.count")
plotDensity(mzXML,z="precursor.mz")

mgf <- readMgfData("szwk901104i19801xms1.mgf")
extractSpectraData(mgf)
methods(class="MSpectra")
mz(mgf)
intensity(mgf)
rtime(mgf)
precursorMz(mgf)
precursorCharge(mgf)
precScanNum(mgf)
precursorIntensity(mgf)
acquisitionNum(mgf)
scanIndex(mgf)
peaksCount(mgf)
msLevel(mgf)
tic(mgf)
ionCount(mgf)
collisionEnergy(mgf)
fromFile(mgf)
polarity(mgf)
smoothed(mgf)
centroided(mgf)
isCentroided(mgf)
writeMgfData(mgf, con = "spectra.mgf", COM = NULL, TITLE = NULL)
removePeaks(mgf, t, msLevel., ...)
filterMsLevel(mgf, msLevel=2)
as.ExpressionSet(mgf)

# This turned to be really slow!
sp_list <- lapply(seq_along(mgf), function(i) {
  intensity_i <- intensity(mgf)[[i]]
  mz_i <- mz(mgf)[[i]]
  centroided_i <- centroided(mgf)[[i]]
  return(new("Spectrum1", intensity = intensity_i, mz = mz_i, centroided = centroided_i))
})
sp1 <- do.call(rbind, sp_list)
# only the first one is more manageable
sp1 <- new("Spectrum1",intensity=intensity(mgf)[[1]],mz=mz(mgf)[[1]],centroided=centroided(mgf)[[1]])
sp2 <- pickPeaks(sp1)
intensity(sp2)
plot(mz(sp1),intensity(sp1),type="h")
## Without m/z refinement
points(mz(sp2), intensity(sp2), col = "darkgrey")
## Using k = 1, closest signals
sp3 <- pickPeaks(sp1, refineMz = "kNeighbors", k = 1)
points(mz(sp3), intensity(sp3), col = "green", type = "h")
## Using descendPeak requiring at least 50% or the centroid's intensity
sp4 <- pickPeaks(sp1, refineMz = "descendPeak", signalPercentage = 50)
points(mz(sp4), intensity(sp4), col = "red", type = "h")


legacy <- function()
{
# individually
  for (raw_file in raw_files) {
      cat("Processing:", raw_file, "\n")
      raw_index <- rawrr::readIndex(raw_file)
      for (scan in raw_index$scan) {
          cat("Scan Number:", scan, "\n")
          raw_spectrum <- rawrr::readSpectrum(raw_file, scan = scan)
          print(raw_spectrum)
      }
  }
  plot(raw_spectrum[[1]], centroid=FALSE)
# a single mzML file
  mzML_file <- list.files(spectra, pattern="mzML",full.names=TRUE, recursive = TRUE)
  dirname(mzML_file)
  mz <- openMSfile(mzML_file)
  fileName(mz)
  runInfo(mz)
  close(mz)
# netCDF, mzXML or mzML
}
