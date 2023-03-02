# Functional activities in daily life

Calcualting functional activities of the upper limb in real life is a challange. Recently Lum et al., (2022) developed a machine learning pipeline to seperate functional from non-functional arm activities in healthy participants and on the non affected side of stroke patients. We tested this model and found that, while it was not yet perfecect, it was better than the previously used activity counts method. 
Therefore, this code will implement the model developed by Lum et al., and tested by Vets et al, on women pre- and post-breast cancer treatment to gain insigt into their arm use in daily life.

### Steps to take before you run this code

1. have matlab installed
2. have the respository downloaded
3. Configure Your System to Use Python
4. This includes setting python to environment variables in windows system path

>  **Note** 
> To call Python® modules in MATLAB®, you must have a supported version of the reference implementation (CPython) installed on your system. 
> This code is created using MATLAB 2022b and Python 3.9 in [Anaconda Distribution](https://www.anaconda.com/products/distribution)
> for other versions of either MATLAB or Python check the [version compatibility](https://nl.mathworks.com/support/requirements/python-compatibility.html)

4. Have the python modules installed

### What does the code do? 

**Data segmentating steps:**
1. Read the native actigraph data usting the pygt3x module in python (implemented in MATLAB).
2. Calculate wear/non-wear time  based on the hip accelerometer.
    - Code from [Syed et al., 2020](https://www.nature.com/articles/s41598-020-62821-2), who implemented [van Hees' Algorithm 1993](https://pubmed.ncbi.nlm.nih.gov/21829556/) using raw acceleration data. We used a non-wear time of 135 minutes and a sliding window of 1 minute, as shown by [Syed et al (2022)](https://www.nature.com/articles/s41598-021-87757-z) to have a high f1-score without needing to resort to deep learning. 
3. If there is at least 5 days and 12 hour per day/block of wear time, than the pre-processing steps start. 

**pre-processing steps:**

1. Refedine axis definition to match those of [Lum et al., (2020)](https://journals.sagepub.com/doi/full/10.1177/1545968320962483)
2. resample data from 30Hz to 50Hz
3. calulate features


**processing**
4. use the pretrained model to predict functional (label 1) and non-functional (label 0) activities of the upper limb
5. calculate the minutes functional active per wear block. 
6. Calulcate the percentage of minutes functional active with respect to the total wear time
7. export all outcome to excel


