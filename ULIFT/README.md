<h1> Kinematical analysis of the ULIFT task </h1>

<h2> Project abstract </h2>

The ULIFT task consists of 4 distinct phases: Phase 1) is from the upper level to the middle level, Phase 2) from the middle level to the floor, Phase 3) from the floor to the middle level, and Phase 4) from the middle level to the upper level. 
These different tasks require a different neuromuscular control in terms of the displacement of the weights—where Phases 1 and 2 can be considered eccentric and Phases 3 and 4 concentric activities. Since the movement of the weight to and from the floor require activities from the trunk and legs alongside upper extremity activity these will not be taken into consideration for the rest of this analysis (too much variability possible, and no sensor information from the lower extremities).

A rough distinction of the Phases 1 and 4 from Phases 2 and 3 happen through the abrupt changes in average position data from the lower arm (Figure 1 upper plot). More precise segmentation of Phases 1 and 4 happen using the sensor acceleration and angular velocity of the lower arm sensor. Phase 1 starts when the rate of change of the Euclidian norm of the lower arm angular velocity is at its first maximum. Similarly, Phase 4 ends with the last prominent minima of the rate of change of the Euclidian norm of the lower arm angular velocity (Figure 1, middle plot).

Detection of the end of phase 1 and start of phase 4 happen through local minima and maxima of the sensor acceleration in X direction. The most prominent local minima correspond to the initiation to move the weight. The less prominent local minima correspond to the release of the weight. Using the knowledge that the height of the hand sensor changes, we select that second to last local minima prior to the substantial change of arm height as the end of phase 1—which corresponds to a less prominent minima, the last prominent minima equal the initiation to phase 2. Similarly, the first most prominent local minima of the sensor’s acceleration signal after the arm's position is high again is selected as the start of phase 4 (Figure 1, lower plot). 

After segmenting the first and fourth phase, we found that there was too much inter-person variability for time curve comparison. Therefore, we only selected the middle arm movement—i.e., grabbing of the second weight to grabbing the subsequent weight—for further processing. This further segmentation was again using the acceleration in X direction.

![(/ULIFT_BC/ULIFT/Uitleg/As expected.jpg)](https://github.com/jillemmerzaal/ULIFT_BC/blob/Version-2/ULIFT/Uitleg/As%20expected.jpg)


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

<h3> Manual plot input </h3>

If the data is as expected, than the code will run automatically. It will determine the start and end points of the seperate phases of the ULIFT task and determine the middle movement of phase 1 and phase 4 as desribed in the project description
