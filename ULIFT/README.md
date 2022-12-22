<h1> Kinematical analysis of the ULIFT task </h1>

<h2> Project abstract </h2>

<h2> How to use the code </h2>

<h3> Instalation </h3>

Download the code from this GitHub page and put in in a folder that is ealiy accesible from your matlab directory. 

> **Warning** 
> Keep the file structure as is

<h3> Workflow </h3>

To run the code, the folder where the code is should be set as the current directory. 
Easiest way to do that is to `cd` in the commant window: 

```
cd C:\Users\u0117545\Documents\GitHub\ULIFT_BC\ULIFT 
```

If you set the correct directory, simply open the main file `ULIFT_preprocessing_V2.m` with the following code in the command window:

```
open ULIFT_preprocessing_V2.m
```

Within that file, you'll need to set the path to the `code` (i.e. current directory) and the `data` on line 44 and 45 respectively. 

The `timepoint` you'd want to analyse on line 56. 
Lastly, indicate if you want to `plot` en check the data where `0 == NO` vs `1 == YES`

```
path.code   = "C:\Users\u0117545\Documents\GitHub\ULIFT_BC\ULIFT";
path.root   = 'C:\Users\u0117545\KU Leuven\An De Groef - DATA';

Timepoint   = 'T0';
movement    = "ULIFT";
plot_or_not = 1;
```

As a final step, give the range of subjects you want to run. 

```
subj = (1:50)
```

After setting all this up, the code should run semi-automatically for the UPLIFT-BC project.
