# X-Ray Absorption Spectroscopy (XAS.py)

## Installation

- Download the 'XAS.py' file and add it to your path folder.
- Include the file in your .py file by using 'from . import XAS as xas'

## Use

- The ICD4 object exports data collected at beamline 4-ID-C at the Advanced Photon Source at Argonne National Laboratory.
- Syntax for use is as follows:
```python
from . import XAS as xas

dire = "./data/4-ID-C/my_data/" #The directory that contains your data
base = "MyData" 				#The base name for all your files (e.g. MyData.0001)

my_sample = xas.IDC4(dire, base, start="1", end="2", shortname = "My Sample")
```
This will load filed 'MyData.0001' and 'MyData.0002' to the IDC4 object and process the data into my_sample.processed_dataframe. 'shortname' defines the legend label for this sample.

To plot data:
```python
fig = plt.figure(1, figsize=(16, 4))

my_sample.plot('STD', color='black')
plt.xlabel(r'Energy / eV', fontsize=15)
```
This plots the 'STD' or Standard data. Other options include 'TEY' (Electron Yield), 'TFY' (Fluorescence Yield), and 'sTFY' (Smoothed Fluorescence Yield).

## Example
### Loading each sample edge into a IDC4 object

```python
dire = "./data/4-ID-C/my_data/"
base = "MyData" 
sample_a = xas.IDC4(dire, base, start="248", end="250", shortname ='Sample A') 
sample_b = xas.IDC4(dire, base, start="244", end="246", shortname = 'Sample B')
```
### Plotting STD data
```python
fig = plt.figure(1, figsize=(16, 4))

sample_a.plot('STD', color='black')
sample_b.plot('STD', color='red')
plt.xlabel(r'Energy / eV', fontsize=15)

show()
```
![Before Alignment](./examples/images/BeforeAlign.jpg "Before Alignment")
### Aligning STD data to account for beamdrift

```python
#Find the peak max value between 705 eV and 710 eV and set the x value to 'x_target'
x_target, y = sample_a.max_in_range('STD', low=705, high=710, plot=False, do_return=True)

#Assign the sample_b peak max value to x1 within the same range
x1, y1 = sample_b.max_in_range('STD', low=705, high=710, plot=False, do_return=True)

#Align the sample_b object to sample_a
sample_b.align(x1, x_target)

#Replot STD data
fig = plt.figure(1, figsize=(16, 4))

sample_a.plot('STD', color='black')
sample_b.plot('STD', color='red')
plt.xlabel(r'Energy / eV', fontsize=15)

show()
```
![After Alignment](./examples/images/AfterAlign.jpg "After Alignment")
