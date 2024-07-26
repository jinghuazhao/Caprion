import numpy as np
from scipy.signal import convolve
from pyteomics import mzml
import csv

def ms_boxcar(mzml_file, charge_states, mass_tolerance, num_peaks=3, min_intensity=0.1):
    """
    MS-BoxCar algorithm for peptide identification in MS/MS spectra from mzML files.

    Parameters:
    - mzml_file (str): Path to mzML file
    - charge_states (list of int): Expected charge states (e.g., [2, 3])
    - mass_tolerance (float): Tolerance in Da (e.g., 0.5)
    - num_peaks (int, optional): Number of peaks to consider (default: 3)
    - min_intensity (float, optional): Minimum intensity threshold (default: 0.1)

    Yields:
    - spectrum_id (str): Scan ID
    - retention_time (float): Retention time (seconds)
    - scores (array): BoxCar scores for each peak
    - peak_indices (array): Indices of top-scoring peaks
    - precursor_mz (float): Precursor m/z
    - precursor_charge (int): Precursor charge
    """
    with mzml.read(mzml_file) as reader:
        for spec in reader:
            if spec['ms level'] == 2:  # Check for MS/MS spectra
                try:
                    # Access precursor information from 'precursorList'
                    precursor_mz = spec['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z']
                    precursor_charge = spec['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['charge state']

                    # Access m/z and intensity arrays
                    mz = spec['m/z array']
                    intensities = spec['intensity array']

                    # Preprocess spectrum
                    if intensities.max() > 0:
                        intensities = intensities / intensities.max()  # Normalize
                    else:
                        print(f"Warning: Maximum intensity is zero for spectrum {spec['id']}. Skipping normalization.")

                    mask = intensities >= min_intensity
                    mz, intensities = mz[mask], intensities[mask]

                    # Compute BoxCar scores for each charge state
                    scores = np.empty((len(charge_states), num_peaks))  # Initialize as NumPy arrays
                    peak_indices = np.empty((len(charge_states), num_peaks), dtype=int) 
                    for i, charge in enumerate(charge_states):
                        scores[i], peak_indices[i] = _compute_boxcar_scores(
                            mz, intensities, charge, mass_tolerance, num_peaks
                        )

                    spectrum_id = spec['id']
                    retention_time = spec['scanList']['scan'][0]['scan start time']

                    # Save scored peaks to a file (e.g., in MGF format)
                    with open(f"spectrum_{spectrum_id.replace('=', '_').replace(' ', '_')}.mgf", "w") as f:  # Sanitize filename
                        f.write(f"BEGIN IONS\n")
                        f.write(f"TITLE=Spectrum {spectrum_id}\n")
                        f.write(f"PEPMASS={precursor_mz}\n")
                        f.write(f"CHARGE={precursor_charge}+\n")
                        for mz_val, intensity_val in zip(mz[peak_indices], intensities[peak_indices]):
                            f.write(f"{mz_val} {intensity_val}\n")
                        f.write(f"END IONS\n")

                    yield (
                        spectrum_id,
                        retention_time,
                        scores,
                        peak_indices,
                        precursor_mz,
                        precursor_charge,
                    )
                except (KeyError, IndexError) as e:
                    print(f"Error processing spectrum {spec['id']}: {e}")
                    continue  # Skip to the next spectrum

def _compute_boxcar_scores(mz, intensities, charge, mass_tolerance, num_peaks):
    """Internal function for computing BoxCar scores."""
    print(f"Intensities: {intensities}")
    print(f"Intensities Length: {len(intensities)}")

    if intensities.size > 0 and not np.all(np.isnan(intensities)):  # Check for empty or all NaN array
        # Calculate mass differences only if mz has more than one element
        if len(mz) > 1:
            mass_diff = np.diff(mz)
        else:
            mass_diff = np.array([1.0])  # Set a default value if mz has only one element

        # Create a boxcar filter with an upper limit on width
        boxcar_width = int(mass_tolerance / np.median(mass_diff))
        boxcar_width = min(boxcar_width, len(intensities) - 1)  # Limit boxcar width

        boxcar_filter = np.ones(boxcar_width)

        # Convolve with intensities
        convolved = convolve(intensities, boxcar_filter, mode='same')

        # Find peaks and their indices
        peak_indices = np.argpartition(convolved, -num_peaks)[-num_peaks:]
        peak_indices = peak_indices[np.argsort(convolved[peak_indices])][::-1]  # Sort by score
        scores = convolved[peak_indices]

        return scores, peak_indices
    else:
        print(f"Warning: Empty or all NaN intensities array for spectrum.")  # Handle empty or all NaN array

# Example usage
mzml_file = 'szwk901104i19801xms1.mzML'  # Your mzML file
charge_states = [2, 3]  # Adjust if necessary
mass_tolerance = 0.5  # Adjust if necessary

with open("ms_boxcar_results.csv", "w", newline="") as csvfile:
    csv_writer = csv.writer(csvfile)
    csv_writer.writerow(["Spectrum ID", "Retention Time", "Scores", "Peak Indices", "Precursor m/z", "Precursor Charge"])

    for result in ms_boxcar(mzml_file, charge_states, mass_tolerance, num_peaks=2, min_intensity=0.2):
        spectrum_id, retention_time, scores, peak_indices, precursor_mz, precursor_charge = result

        # Convert peak_indices to a string representation (now works correctly)
        peak_indices_str = np.array2string(peak_indices, separator=',')

        csv_writer.writerow([spectrum_id, retention_time, scores, peak_indices_str, precursor_mz, precursor_charge])
