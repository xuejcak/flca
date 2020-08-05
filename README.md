# FLCA
Fourier local correlation alignment (FLCA) is an efficient method for aligning images. The core of the code is changed from [FLCT](http://solarmuri.ssl.berkeley.edu/overview/publicdownloads/software.html) (Fisher & Welsch PASP 383, 373, (2008), Fisher et al. ApJS 248, 2, (2020)), and the Python code is written by Xue,J. (xuejc@pmo.ac.cn). The FLCT is a method for performing local correlation tracking. FLCA calculates the local corretion of two images pixel by pixel, and the final shift-X and shift-Y are gained as the median values of those from the central region. The code has the advantages of being precise, fast, user friendly, and open.

## License
The FLCA is open-source, under the lesser GPL version 2.1 license. You are free to copy the code and use it as you like.

## Installation and run
The core of the code is changed from the FLCT, and is written in C, which is called by the main Python code. The library fftw3 and Python modules numpy, astropy, skimage, ctypes are necessary. After editting the ./source/Makefile for the fftw3 path, you may enter the ./source directory and compile the C file in the terminal:
```bash
cd (where flca.py is)/source
make
```
After that, there should be a new file ./source/flca.so.
For running the code, you and use the GUI, or edit the input.txt and run flca using Python3. For the first mode, run the following code in the terminal
```bash
python3 flca_gui.py
```
or run flca_gui in Python IDLE, Jupyter notebook, spyder, etc.. You need to enter the required information and click "Single-channel coalignment" or "Double-channel coalignment". 
For the second mode, you need to enter the information in the "input.txt" (or using other filename), and run the following codes in Python.
```python
import os
os.chdir('The directory where flca.py is')
import flca
flca.run('./input.txt') # or other filename.
```
If you coalign files to those in the same channel one by one, the coalignment error for the nearby two files are small, but the files tend to move to a certain direction after coaligning tens of files. Then you need to sample some files and coalign them with other filtergrams, such as images from SDO/AIA. Record the file index (ref) and shifts (xm and ym) in detrend.py, and run it:
```bash
python3 detrend.py
```
or
```python
detrend
```
 The sampled files will be shifted as defined xm and ym, and other files are shifted linearly.

## Parameters and keywords
FLCA read the necessary information throught the GUI or "input.txt". The FLCA offer two choices for you, one is that files are coaligned to those in the same channel (Default, "Single-channel coalignment" in GUI), and the other is that files are coaligned to those beloning to the other channel (twochannel in "input.txt" or "Double-channel coalignment" in GUI). In the latter case, the files in the two channels must have the same pixel scale.
In the single-channel case, the following parameters and keywords are available.
- indir: The directory that contains the files you want to align. 
- outdir: The directory that the aligned files will be written to.
- firstfits: The first reference file, which must have the same name as one of the fits files in indir. Different directories are allowed. Default is the first fits file in the "indir".
- x0, x1, y0, y1: Indices of the pixels to define the region that you want to compare.
- every: After "every" number files, the reference file will be changed. Default is 1. Even "every" is set to be >1, I find that the files still move to one direction (for two limb cases). In this case, every=1 is recommended and detrend.py should be used.
- skip: Generally one does not need to calculate shifts at every location (pixel), set "skip" means that every "skip" pixels in each direction will calculate once shifts.
- xoffset, yoffset: If skip is set, the calculating starts at xoffset and yoffset. Default 0, 0.
- sigma: The standard deviation of the Gaussian mesh, which also determines the size of the mesh. sigma value must be larger than the maximum shift desired, but the larger the sigma is, the longer the running time will be.
- threshold: If the mean value of the two images is lower than that value, the calculation at this location will be cancelled. If the value ends with "a", then it means absolute value. Otherwise, the value between 0 and 1 will be treated as a relative value.
- kr: The maximum wavenumber of the low pass filtering. 0.2-0.5 is recommended.
- biascorrect: The calculated shifts tend to be smaller. If this keyword is set, bias correciton is performed.
- interpolate: If set, shifts are linearly interpolated if shift-X and shift-Y are decimals. Otherwise, the nearest values are adopted.

In the double-channel mode, "firstfits" and "every" are ignored. Meanwhile, the following parameters are necessary.
- indir0: The directory containing the reference files.
- time0a, time0b: The string indices from the beginning of the date (time0a) to the end of the time (time0b) in the filename of references. The negative indices are recommended. For example, if filename[-23:-8] contains the datetime, then set time0a to be -23 and time0b to be -8.
- time1a, time1b: The string indices for the processing files. If not set, they would be the same as time0a and time0b.
