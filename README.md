# HeartLab Data Analysis Workflow

**Author:** Tainan Cerqueira Neves - HEartLab

## General Information

1. Never work directly in the folder if you have a local GIT repository.
   - This could generate modifications in the main branch and affect all repository users.
   - Place the files to be used in your current working folder and make modifications there.
   - As it is a collaborative repository, changes and improvements to the user experience can be made and committed for everyone to access updates.

2. This step-by-step guide is developed to assist collaborative data processing within the HeartLab group, ensuring continuous and facilitated learning regarding signal processing activities within the group.

3. The "00 – examples" folder contains example files used for learning how to use the codes. Due to their size, they should be downloaded to a cloud folder available at: https://1drv.ms/f/s!AmRGmGB_stE1iddm6pYJEiSEWqlHHQ?e=gihHmW

## 01 – seq_to_mat: Binning and Optical Data Extraction

For optical data manipulation, significant processing power is required. Since most available machines lack this power, we perform a procedure called Binning, reducing the size and pixel quantity of each frame, making optical files lighter.

1. Convert .SEQ files (Norpix) to .MAT files (Matlab).
2. Plot and save demonstrative figures of the camera and signal at a point for each recording.
3. Binning – Aggregates pixels into a single pixel by averaging and transforms the color scale to 16 bits.
4. Export binned data.

### How to Use:

1. Place the codes within your working folder.
2. Original optical data should be in a folder named "Optical Mapping" -> "Camera 1" + "Camera 2" + "Camera 3".
3. Create a folder named "Optical Mapping to MATLAB".
4. In MATLAB, right-click on the "Optical Mapping" folder -> "Add to Path" -> "Selected Folders and Subfolders".
   - Do this for the folders "Optical Mapping," "Common Matlab Function," and "Optical Mapping to Matlab".
5. Click to run the code.

**Note:**
- This code requires substantial processing capacity and takes hours to complete. It is recommended to run it overnight. Monitor the first exported files to ensure the code is running correctly.

## 02 – extraction_filtering: Data Extraction and Filtering

These codes perform the extraction of electrical data and the filtering of electrical and optical data.

### How to Use:

1. Optical filtering (..\02 - extraction_filtering\optic):
   - Run "main_optical_filtering.mat".
   - Note: Some code segments should be run individually for all cameras before exporting.
   - The code is well-documented.
   - Use Plot implay to visualize and assess recording quality.
   - Generates a .MAT file with filtered optical data and corresponding ROIs.

2. Electrical filtering (..\02 - extraction_filtering\electric\ open_ephys_extraction_matlab_filtering):
   - Run "main_extraction_filtering.mat".
   - The code is self-explanatory; follow the step-by-step instructions.
   - This code allows immediate data visualization, so use it to analyze your data quality.
   - Exports two files: one with raw electrical data and another with filtered electrical data.

## 03 - synchronization_optical_electric: Synchronization of Optical and Electrical Data

Here, optical data is resampled to a 4000 Hz frequency to match the electrical frequency, facilitating future analyses and plots.

### How to Use:

1. Run the code "main_synchronization" (..\03 - synchronization_optical_electric).
   - The code exports a variable containing all synchronized optical and electrical data.

## 04 – signal_plot: Signal Plots

### How to Use:

1. Run the code "main_signal_plot.mat" (..\04 - signal_plot).
   - It contains various plot options.

## 05 – potential_maps: Potential Maps

There are three folders within this directory, two related to maps for electrical data and one for optical data.

### How to Use:

1. Enter the folder for the code you want to use and run "main_".
   - For electrical – 3D_mode – Laplacian: "main_electric_potential_plot.m"
   - For electrical – matrix_mode – close neighbors: "main"
   - For optic: "main_optic_potential_map.m"

## 06 - dominant_frequency_maps: Dominant Frequency Maps

### How to Use:

1. Run the code "main_dominant_frquency_map.m" (..\06 – dominant_frequency_maps).
   - Run it point by point using F9. It is well-documented and generates dominant frequency plots for MEAs, TANK, and OPTIC.

## 07 – phase_maps: Phase Maps

### How to Use:

1. Run the code "main_phase_map.m" (..\07 – phase_maps).
   - Run point by point using F9 to ensure control of occurrences.

## 08 – 

# To Do:

1. Add lines demarcating the faces of the tank in plots containing a tank.
2. DONE - Change electrode indicators to black (currently in red).
3. In the optical export code, create a figure that is load with low filter to guide electrodes position. Keep the export file smooth, but export a single frame to guide this selection.
4. DONE - Simplify the Potential plot code. Use the modularization method
5. DONE - In potential plot, add the layer with the vertices
6. Create a switch to select cameras for filtering instead of set it manually
7. Verify why there is one more sample after syncronization and correct it
8. Save the pixel position by clicking instead of insert it manually (Pick up a trace) - Save in variable
9. Add CL code
for i=1:(size(A,1)-1)
   s=A(i,:);
   [y,x]=findpeaks(s,'MinPeakHeight',0.005,'MinPeakDistance',1000);
   if ~isempty(x)
      cl=x(2)-x(1);
      CL(i)=cl;
   end
end