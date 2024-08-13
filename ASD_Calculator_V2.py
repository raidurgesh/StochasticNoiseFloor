'''
Owen Coffee
8/12/24
ASD_Calculator_V2.py
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import welch, find_peaks
import os

# Parameters
Time = 10000
Ntime = Time + 1

# Generate random noise data (Uniform distribution)
time_s = np.linspace(0, Ntime, num=Ntime, endpoint=False)
array = np.random.rand(Ntime)
rand_p1K = 0.1 * (array * 2. - 1.)
data_p1K = 20. + rand_p1K

# Generate Gaussian (Normal) noise data
mean = 0
stdp1 = 0.04
Normalp1 = np.random.normal(mean, stdp1, size=Ntime)
Ndata_p1K = 20. + Normalp1

# Sort the uniform and Gaussian noise data
ABp1K = np.abs(data_p1K)
SHp1K = np.sort(ABp1K)

AB_Ndata_p1K = np.abs(Ndata_p1K)
SH_Ndata_p1K = np.sort(AB_Ndata_p1K)

# Compute the Power Spectral Density (PSD) using Welch's method
fs = 1.0  # Sampling frequency (Hz), assuming 1 Hz for simplicity
frequencies, psd_uniform = welch(data_p1K, fs=fs, nperseg=1024)
frequencies, psd_sorted_uniform = welch(SHp1K, fs=fs, nperseg=1024)

frequencies, psd_gaussian = welch(Ndata_p1K, fs=fs, nperseg=1024)
frequencies, psd_sorted_gaussian = welch(SH_Ndata_p1K, fs=fs, nperseg=1024)

# Calculate the Amplitude Spectral Density (ASD) from the PSD
asd_uniform = np.sqrt(psd_uniform)
asd_sorted_uniform = np.sqrt(psd_sorted_uniform)

asd_gaussian = np.sqrt(psd_gaussian)
asd_sorted_gaussian = np.sqrt(psd_sorted_gaussian)

# # Find peaks in the PSD and ASD
# peaks_psd_uniform, _ = find_peaks(psd_uniform)
# peaks_asd_uniform, _ = find_peaks(asd_uniform)
# peaks_psd_sorted_uniform, _ = find_peaks(psd_sorted_uniform)
# peaks_asd_sorted_uniform, _ = find_peaks(asd_sorted_uniform)

# peaks_psd_gaussian, _ = find_peaks(psd_gaussian)
# peaks_asd_gaussian, _ = find_peaks(asd_gaussian)
# peaks_psd_sorted_gaussian, _ = find_peaks(psd_sorted_gaussian)
# peaks_asd_sorted_gaussian, _ = find_peaks(asd_sorted_gaussian)

# Plot the results
plt.figure(figsize=(18, 21))

# Plot the uniform noise data
plt.subplot(7, 1, 1)
plt.plot(time_s, data_p1K, label='Uniform Noise')
plt.title('Uniform Random Noise Data')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.grid(True)
plt.legend()

# Plot the Gaussian noise data
plt.subplot(7, 1, 2)
plt.plot(time_s, Ndata_p1K, label='Gaussian Noise', color='orange')
plt.title('Gaussian Noise Data')
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')
plt.grid(True)
plt.legend()

# Plot the PSD comparison
plt.subplot(7, 1, 3)
plt.loglog(frequencies, psd_uniform, label='Uniform Noise PSD')
plt.loglog(frequencies, psd_gaussian, label='Gaussian Noise PSD', color='orange')
# plt.plot(frequencies[peaks_psd_uniform], psd_uniform[peaks_psd_uniform], 'rx', label='Uniform PSD Peaks')
# plt.plot(frequencies[peaks_psd_gaussian], psd_gaussian[peaks_psd_gaussian], 'bx', label='Gaussian PSD Peaks')
plt.title('Power Spectral Density (PSD) Comparison')
plt.xlabel('Frequency (Hz)')
plt.ylabel('PSD (W/Hz)')
plt.grid(True, which='both', ls='--')
plt.legend()

# Plot the ASD comparison
plt.subplot(7, 1, 4)
plt.loglog(frequencies, asd_uniform, label='Uniform Noise ASD')
plt.loglog(frequencies, asd_gaussian, label='Gaussian Noise ASD', color='orange')
# plt.plot(frequencies[peaks_asd_uniform], asd_uniform[peaks_asd_uniform], 'rx', label='Uniform ASD Peaks')
# plt.plot(frequencies[peaks_asd_gaussian], asd_gaussian[peaks_asd_gaussian], 'bx', label='Gaussian ASD Peaks')
plt.title('Amplitude Spectral Density (ASD) Comparison')
plt.xlabel('Frequency (Hz)')
plt.ylabel('ASD (W/√Hz)')
plt.grid(True, which='both', ls='--')
plt.legend()

# Plot the ASD comparison (Uniform Noise)
plt.subplot(7, 1, 5)
plt.loglog(frequencies, asd_uniform, label='Original Uniform Noise ASD')
plt.loglog(frequencies, asd_sorted_uniform, label='Sorted Uniform Noise ASD', color='green')
plt.title('Amplitude Spectral Density (ASD) - Uniform Noise')
plt.xlabel('Frequency (Hz)')
plt.ylabel('ASD (W/√Hz)')
plt.grid(True, which='both', ls='--')
plt.legend()

# Plot the ASD comparison (Gaussian Noise)
plt.subplot(7, 1, 6)
plt.loglog(frequencies, asd_gaussian, label='Original Gaussian Noise ASD')
plt.loglog(frequencies, asd_sorted_gaussian, label='Sorted Gaussian Noise ASD', color='red')
plt.title('Amplitude Spectral Density (ASD) - Gaussian Noise')
plt.xlabel('Frequency (Hz)')
plt.ylabel('ASD (W/√Hz)')
plt.grid(True, which='both', ls='--')
plt.legend()

# Plot the sorted noise data for comparison
plt.subplot(7, 1, 7)
plt.plot(SHp1K, label='Sorted Uniform Noise', color='green')
plt.plot(SH_Ndata_p1K, label='Sorted Gaussian Noise', color='red')
plt.title('Sorted Noise Data')
plt.xlabel('Sample Index')
plt.ylabel('Amplitude')
plt.grid(True)
plt.legend()

plt.tight_layout()

# # Ask the user if they want to save the file
# save_file = input("Do you want to save the plots and peak data to a PDF file? (yes/no): ").strip().lower()

# if save_file == 'yes':
#     file_path = input("Please enter the file path where you want to save the PDF (including the file name, e.g., 'output.pdf'): ").strip()
    
#     # Save the figure as a PDF
#     plt.savefig(file_path)

#     # Also, save peak information in a text file
#     base_path, _ = os.path.splitext(file_path)
#     peak_info_path = base_path + "_peaks.txt"
    
#     with open(peak_info_path, 'w') as f:
#         f.write("Peaks in Original Uniform Noise PSD:\n")
#         for i in peaks_psd_uniform:
#             f.write(f"Frequency: {frequencies[i]} Hz, PSD: {psd_uniform[i]} W/Hz\n")
        
#         f.write("\nPeaks in Sorted Uniform Noise PSD:\n")
#         for i in peaks_psd_sorted_uniform:
#             f.write(f"Frequency: {frequencies[i]} Hz, PSD: {psd_sorted_uniform[i]} W/Hz\n")
        
#         f.write("\nPeaks in Original Gaussian Noise PSD:\n")
#         for i in peaks_psd_gaussian:
#             f.write(f"Frequency: {frequencies[i]} Hz, PSD: {psd_gaussian[i]} W/Hz\n")
        
#         f.write("\nPeaks in Sorted Gaussian Noise PSD:\n")
#         for i in peaks_psd_sorted_gaussian:
#             f.write(f"Frequency: {frequencies[i]} Hz, PSD: {psd_sorted_gaussian[i]} W/Hz\n")
        
#         f.write("\nPeaks in Original Uniform Noise ASD:\n")
#         for i in peaks_asd_uniform:
#             f.write(f"Frequency: {frequencies[i]} Hz, ASD: {asd_uniform[i]} W/√Hz\n")
        
#         f.write("\nPeaks in Sorted Uniform Noise ASD:\n")
#         for i in peaks_asd_sorted_uniform:
#             f.write(f"Frequency: {frequencies[i]} Hz, ASD: {asd_sorted_uniform[i]} W/√Hz\n")
        
#         f.write("\nPeaks in Original Gaussian Noise ASD:\n")
#         for i in peaks_asd_gaussian:
#             f.write(f"Frequency: {frequencies[i]} Hz, ASD: {asd_gaussian[i]} W/√Hz\n")
        
#         f.write("\nPeaks in Sorted Gaussian Noise ASD:\n")
#         for i in peaks_asd_sorted_gaussian:
#             f.write(f"Frequency: {frequencies[i]} Hz, ASD: {asd_sorted_gaussian[i]} W/√Hz\n")
    
#     print(f"Plots saved to {file_path}")
#     print(f"Peak information saved to {peak_info_path}")

plt.show()
