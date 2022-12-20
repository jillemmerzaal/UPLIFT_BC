<h1 align="center">Movement Quality</h1>

<p align="center">
Code calculate movement quality to determine upper limb function in people before and after breast cancer surgery
</p>

<h2> Project description </h2>
This project is written with the single purpose of working for the UPLIFT-BC project.

> **Note**
> For any other function or purpose, this code needs to be addapted to fit. 

Within this project, the following movement quality parameters will be calculated per participant, per activit: 
1) dynamic stability of movement as defined by the Lyapunov Exponent
2) signal smoothness as defined by the log dimensionless jerk
3) signal predictability with sample entropy
4) signal regularity with the autocorrelation
5) movement speed based on the repetitions
6) movement variability based on the repetions
7) movement variabilty from root mean square and root mean square ratio

The calculated movement quality parameters will be stored in an Excel file in the Output folder within this repository. IF the output file excists, the data will be appended. ELSE a new file will be created. 
> **Warning** 
> Keep in mind that the data will be added. Thus if you run a subject twice, the data will be stored twice!


<h2> Instalation </h2>
Download the code from this GitHub page and put in in a folder that is ealiy accesible from your matlab directory. 

> **Warning** 
> Do NOT move the files into seperate folders, and keep the structure as is. 


To run the code, the folder where the code is should be set as the current directory. 
Easiest way to do that is to **cd <path to your directory>** in the commant window: 

```
cd C:\Users\u0117545\Documents\GitHub\ULIFT_BC\MovQual 
```

<h2> How to run </h2>
