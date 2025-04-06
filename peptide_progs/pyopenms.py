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

import pyopenms as oms
import numpy as np
import pandas as pd
import datashader as ds
import holoviews as hv
import holoviews.operation.datashader as hd
from holoviews.plotting.util import process_cmap
from holoviews import opts
import sys

exp = oms.MSExperiment()
loader = oms.MzMLFile()
loadopts = loader.getOptions()
loadopts.setSkipXMLChecks(True)
loadopts.setIntensity32Bit(True)
loadopts.setIntensityRange(oms.DRange1(oms.DPosition1(5000), oms.DPosition1(sys.maxsize)))
loader.setOptions(loadopts)
mzml_file_path = "/home/jhz22/Caprion/analysis/UHZ/mzML/szwk021728f22101xms1.mzML"
loader.load(mzml_file_path, exp)
exp.updateRanges()
ms_level = 1
min_rt = 0.0
max_rt = 1000.0
min_mz = 100.0
max_mz = 2000.0
rt, mz, intensity = exp.get2DPeakDataLong(min_rt, max_rt, min_mz, max_mz, ms_level)
spectra_df = pd.DataFrame({
    'RT': rt,
    'm/z': mz,
    'Intensity': intensity
})
spectra_df.set_index(['RT', 'm/z'], inplace=True)
points = hv.Points(spectra_df, kdims=["RT", "m/z"], vdims=["Intensity"], label="MS1 survey scans").opts(
    fontsize={"title": 16, "labels": 14, "xticks": 6, "yticks": 12},
    colorbar=True,
    width=1000,
    height=1000,
    tools=["hover"]
)
raster = hd.rasterize(points, cmap=process_cmap("blues", provider="bokeh"), aggregator=ds.sum("Intensity"), cnorm="log", alpha=10, min_alpha=0).opts(
    active_tools=["box_zoom"],
    tools=["hover"],
    plot=dict(
        width=800,
        height=800,
        xlabel="Retention time (s)",
        ylabel="mass/charge (Da)"
    )
)
spread_raster = hd.dynspread(raster, threshold=0.7, how="add", shape="square")
spread_raster
