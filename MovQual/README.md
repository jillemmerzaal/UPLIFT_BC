<h1> Movement Quality</h1>

Code to calculate movement quality to determine upper limb function in people before and after breast cancer surgery.

> **Note**
> This project is written with the sole purpose of working for the UPLIFT-BC project. 

<h2> Project abstract </h2>

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

Detailed information about the movement quality parameters can be found at the end of this document.

The calculated movement quality parameters will be stored in an Excel file in the Output folder within this repository. If the output file excists, the data will be appended. Else a new file will be created. 
> **Note** 
> Keep in mind that the data will be added to an existing document. Thus if you run a subject twice, that subject will be stored twice!

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


<h2> Parameters </h2>

<h4> Lyapunov Exponent </h4>

The Lyapunov exponent is a mathematical concept that helps to quantify the rate at which the system becomes regular over time. It is commonly used to study the behaviour of chaotic systems over time, such as weather reports or the stock market. 
In biomechanics, the LyE can be used to study the stability of human movement. For example, when a person is walking or running, small perturbations or disturbances in their movement can cause the center of mass to move in unpredictable ways. By calculating the LyE for these movements, we can determine how stable the movement is and the sensitivity to initial conditions. By understanding how sensitive the movement is to changes in the initial condition, we could develop strategies to improve performance and reduce risk of injury. 
In terms of arm movements, researchers have used the Lyapunov Exponent to study the stability of reaching movements, playing an instrument, or throwing movements. Overall, the Lyapunov exponent is a usefull tool for understangin the behaviour of movement in a variety of settings including sports, rehabilitation and ergonomics.  

<h4> Sample Entropy </h4>

Sample Entropy is a measure of randomness or predictability of a time series data set. It is a mathematical algorithm tha compares the similarity of patterns withing a dataset over different times scales and can be used to assess predictability of a signal. In simple terms, sample entropy lookts at how similar patterns are within a dataset. A data set with a high sample entropy would indicate a high degree of randomness *(unpredictable)*, while a data set with low sample entropy would indicate low degree of randomness *(predictable)*. 
In terms of arm movement, Sample Entropy has been used to assess reaching tasks in people with parkinson's disease compared to control subjects. Results showed that people with parkinson's had a lower sample entropy values, indicating that arm movement was less complex. Concluding that sample entropy could be used as a tool to assess arm movement complexity in individuals with parkinsons's disease and for evaluating the effect of different treatments on arm movement. It has also been used to assess the complexity of arm movements in individuals with cerebellar ataxia and during different stages of development in children. 
In the field of biomechancis, sample entropy have been used to assess the gait or running pattern of an individual to determine if there are abnormalities or changes over time. Where too low and too high sample entropy would indicate a *posible* problem and restrict people in changing their motor control pattern. OVerall, sample entropy is a usefull tool for analysing time series data in the field of biomechanics and can provide valuable insight into verious physiological processess.
