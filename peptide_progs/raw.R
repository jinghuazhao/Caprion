options(width=200)

# ZWK .raw data
library(rawrr)
spectra_ZWK <- "~/Caprion/pre_qc_data/spectral_library_ZWK"
raw_files <- list.files(spectra_ZWK, pattern = "\\.raw$", full.names = TRUE)
## collectively
suppressMessages(library(MsBackendRawFileReader))
ZWKw <- Spectra::backendInitialize(MsBackendRawFileReader::MsBackendRawFileReader(),
        files = raw_files)
class(ZWK)
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
  library(mzR)
  mzML_file <- list.files(spectra, pattern="mzML",full.names=TRUE, recursive = TRUE)
  dirname(mzML_file)
  mz <- openMSfile(mzML_file)
  fileName(mz)
  runInfo(mz)
  close(mz)
# netCDF, mzXML or mzML
}
