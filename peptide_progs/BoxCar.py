import numpy as np
from scipy.signal import convolve

def ms_boxcar(mz, intensities, charge, mass_tolerance, **kwargs):
    """
    MS-BoxCar algorithm for peptide identification in MS/MS spectra.

    Parameters:
    - mz (array): m/z values
    - intensities (array): Corresponding peak intensities
    - charge (int): Charge state (z)
    - mass_tolerance (float): Tolerance in Da
    - **kwargs: Passed to boxcar_kernel (e.g., num_peaks, min_intensity)

    Returns:
    - scores (array): BoxCar scores for each peak
    - peak_indices (array): Indices of top-scoring peaks
    """
    def boxcar_kernel(mass, charge, tolerance, num_peaks, min_intensity):
        """Compute boxcar kernel and apply intensity thresholding"""
        spacing = 1.0 / charge  # Isotopic spacing
        kernel_size = int(2 * tolerance / spacing) + 1
        kernel = np.ones(kernel_size)  # Uniform boxcar
        kernel[0] = kernel[-1] = 0.5  # Half-weight edges
        intensities_thr = intensities >= min_intensity
        return kernel, intensities_thr, kernel_size  # Return kernel_size

    num_peaks = kwargs.get('num_peaks', 3)
    min_intensity = kwargs.get('min_intensity', 0.1)

    kernel, intensities_thr, kernel_size = boxcar_kernel(mz[0], charge, mass_tolerance, num_peaks, min_intensity)
    scores = np.correlate(intensities[intensities_thr], kernel, mode='valid')

    # Identify top-scoring peaks
    idx = np.argsort(-scores)[:num_peaks]
    return scores[idx], idx + (kernel_size - 1) // 2  # Correctly use kernel_size

def identify_peptides(spectra, charges, mass_tolerance, **kwargs):
    """
    Batch processing for multiple MS/MS spectra.

    Parameters:
    - spectra (list of tuples): [(mz, intensities), ...]
    - charges (list of int): Corresponding charge states
    - mass_tolerance (float): Tolerance in Da
    - **kwargs: Passed to ms_boxcar (e.g., num_peaks, min_intensity)

    Returns:
    - results (list of tuples): [(scores, peak_indices), ...]
    """
    results = []
    for (mz, intensities), charge in zip(spectra, charges):
        scores, peak_indices = ms_boxcar(mz, intensities, charge, mass_tolerance, **kwargs)
        results.append((scores, peak_indices))
    return results

spectra = [
    (np.array([100.0, 101.0, 102.0, 103.0]), np.array([10, 20, 30, 40])),
    (np.array([200.0, 201.0, 202.0, 203.0]), np.array([50, 60, 70, 80])),
]
charges = [2, 3]
mass_tolerance = 0.5

results = identify_peptides(spectra, charges, mass_tolerance, num_peaks=2, min_intensity=0.2)

for i, ((scores, peak_indices), (mz, _)) in enumerate(zip(results, spectra)):
    print(f"Spectrum {i+1}:")
    print(f"  Top scores: {scores}")
    print(f"  Corresponding peaks: {mz[peak_indices]}")
    print()
