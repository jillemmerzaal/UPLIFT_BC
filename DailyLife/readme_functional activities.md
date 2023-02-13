
# Functional activities in daily life

Calcualting functional activities of the upper limb in real life is a challange. Recently Lum et al., (2022) developed a machine learning pipeline to seperate functional from non-functional arm activities in healthy participants and on the non affected side of stroke patients. We tested this model and found that, while it was not yet perfecect, it was better than the previously used activity counts method. 
Therefore, this code will implement the model developed by Lum et al., and tested by Vets et al,m on women pre- and post-breast cancer treatment to gain insigt into their arm use in daily life.

### Steps to take before you run this code

1. have matlab installed
2. have the respository downloaded
3. Configure Your System to Use Python

> To call Python® modules in MATLAB®, you must have a supported version of the reference implementation (CPython) installed on your system. 
> This code is created using MATLAB 2022b and Python 3.9 in [Anaconda Distribution](https://www.anaconda.com/products/distribution)
> 
> For other versions of either MATLAB or Python check the [version compatibility](https://nl.mathworks.com/support/requirements/python-compatibility.html)

4. Have the python modules installed || see requirements.txt

### What does the code do? 

**Data segmentating steps:**
1. read the native actigraph data and transport it to csv format
> Use the python code in a matlab environment
2. automatically determine when the signal starts (i.e. watches are put on)
3. calculate wear/non-wear time
4. if the data equivalent of 5 days and 12 hours per day is reached, the process will continue

**pre-processing steps:**



