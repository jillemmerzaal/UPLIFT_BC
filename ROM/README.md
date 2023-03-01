# Kinematical analysis for range of motion 

This is the main repository for the upper limb range of motion timecurve analysis as calculated by the xsens system. Data collected in funtion of the UPLIFT_BC project using Xsens MVN 2021.2 software.

**specify**
- Path where all the code is stored
- Path where all the code sections are stored (if not in path.orig)
- Timepoint you want to analyse (e.g. 'T0' or 'T1')
- movement you want to analyse (i.e. 'Af', 'EXO', or 'ABD')
- path where all the data is stored
- define if you want to check the data: plot_or_not = 1 & check_complete = 0
- define if you want to safe the data: plot_or_not = 0 & check_complete = 1
- define the range of subjects you want to analyse (subj = 1:10 runs the first 10 subjects | subj = 10 only runs subject 10 | subj = (1:5 7:10) runs subject 1 through 10 except 6 ect.)
```
path.code       = "C:\Users\u0117545\Documents\GitHub\ULIFT_BC\ROM";
path.root       = 'C:\Users\u0117545\KU Leuven\An De Groef - DATA';

Timepoint       = 'T1';
movement        = 'ABD';
plot_or_not     = 0;
check_complete  = 1;
scapulo         = 0;
secondary       = 0;

subj            = [2:4, 6:7, 9:10]
```

**Secondary and scapulothoratic angles**
```scapulo = 1```  also extracts the 3d scapulothoratic movement during the primary movement.
```secondary=1``` also include the other two movements in the glenohumeral joint. e.g. if the momvent is ABD (abduction) than primary will only extract the abduction movement from the xsens data. scapulo, will extract abduction from the glenohumeral joint, plus 3D scapulotoracaal angles. secondary will extract the all glenohumaral joint angles. Thus, if you want glenohumeral abduction + the scapulothoratic joint angles than: ```scapulo = 1```; ```secondary = 0```


**What does the code do?**
First, mvnx. files are loaded into matlab, and based on the motion the primary kinimatical timecurve is extracted (i.e. if Motion = 'ABD', the abduction of the shoulder is the primary kinamtical timecurve). 
Than the different repetitions get detected using the findsignal function in matlab. The comparison signal is half a sin wave. 
From there, the seperate repetions are extracted, timenormalised to 101 datapoints (i.e. 100%), averaged over time and written to an excell file. 




