<h1 align="center">Movement Quality</h1>

<p align="center">
Code calculate movement quality to determine upper limb function in people before and after breast cancer surgery
</p>

<h2> Project description </h2>

> **Note**
> This project is written with the single purpose of working for the UPLIFT-BC project. 

Participants performed a weighted pseudo-cyclic motion task. They wore an IMU (Xsens technologies, MVN awinda, 60Hz) on each wrist. The weight they were holding was aproximatly 1 kg. They moved the weight from the side of their body to a shelf that was fixed at shoulder height and back. This was repeated 14 times and performed at two timepoints: baseline and follow-up at 1 month post-surgery. 
From each trial, the raw acceleration signal was extracted, and from that the movement quality parameters were calculated:

Within this project, the following movement quality parameters will be calculated per participant, per activit: 
- dynamic stability || Lyapunov Exponent and state space reconstruction
- signal smoothness || Log Dimensionless Jerk
- signal predictability || Sample Entropy
- signal regularity || Autocorrelation 
- movement variabilty || Root Mean Square and Root Mean Square Ratio
- movement speed 
- movement variability 

The calculated movement quality parameters will be stored in an Excel file in the Output folder within this repository. IF the output file excists, the data will be appended. ELSE a new file will be created. 
> **Warning** 
> Keep in mind that the data will be added. Thus if you run a subject twice, the data will be stored twice!

<h2> Instalation </h2>
Download the code from this GitHub page and put in in a folder that is ealiy accesible from your matlab directory. 

> **Warning** 
> Do NOT move the files into seperate folders, and keep the structure as is. 


To run the code, the folder where the code is should be set as the current directory. 
Easiest way to do that is to `cd` in the commant window: 

```
cd C:\Users\u0117545\Documents\GitHub\ULIFT_BC\MovQual 
```

<h2> How to run </h2>
