from pyopenms import MzMLFile, MSExperiment
exp = MSExperiment()
MzMLFile().load("/home/jhz22/Caprion/analysis/UHZ/mzML/szwk021728f22101xms1.mzML", exp)
ms2_spectra = []
for spectrum in exp.getSpectra():
    if spectrum.getMSLevel() > 1:
        ms2_spectra.append(spectrum)

spectrum = exp.getSpectra()[0]
mz, intensity = spectrum.get_peaks()
import matplotlib.pyplot as plt
plt.figure(figsize=(10, 6))
plt.stem(mz, intensity, basefmt=" ")
plt.xlabel('m/z')
plt.ylabel('Intensity')
plt.title('MS/MS Spectrum')
plt.show()
