<h1> Movement Quality</h1>

Code to calculate movement quality to determine upper limb function in people before and after breast cancer surgery

<h2> Project abstract </h2>

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

The calculated movement quality parameters will be stored in an Excel file in the Output folder within this repository. If the output file excists, the data will be appended. Else a new file will be created. 
> **Note** 
> Keep in mind that the data will be added. Thus if you run a subject twice, the data will be stored twice!

<h2> How to use the code </h2>

<h3> Instalation </h3>

Download the code from this GitHub page and put in in a folder that is ealiy accesible from your matlab directory. 

> **Warning** 
> Keep the file structure as is

<h3> Workflow </h3>

To run the code, the folder where the code is should be set as the current directory. 
Easiest way to do that is to `cd` in the commant window: 

```
cd C:\Users\u0117545\Documents\GitHub\ULIFT_BC\MovQual 
```

If you set the correct directory, simply open the main file `MovQualityF.m` with the following code in the command window:

```
open MovQulaity.m
```

Within that file, you'll need to set the current directory on line 11, the `timepoint` you'd want to analyse on line 13, and the path where the `data` and the `code` is located on line 15 and 16 respectively. 
Lastly, indicate if you want to `plot` en check the data and/or `save_to_excel` where `0 == NO` vs `1 == YES`


```ruby
cd("C:\Users\u0117545\Documents\GitHub\ULIFT_BC\MovQual")

Timepoint       = 'T0';
movement        = "F";
path.data       = 'C:\Users\u0117545\KU Leuven\An De Groef - DATA';
path.root       = 'C:\Users\u0117545\Documents\GitHub\ULIFT_BC\MovQual';
plot_or_not     = 0;
safe_to_excel   = 1;
``` 
As a final step, give the range of subjects you want to run. 

```
subj = (1:50)
```

After setting all this up, the code should run automatically for the UPLIFT-BC project.
