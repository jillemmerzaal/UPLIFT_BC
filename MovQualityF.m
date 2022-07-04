%% movement quality functional activity

clear all; %close all; clc
%% 1. input data
cd("C:\Users\u0117545\Documents\GitHub\ULIFT_BC")
addpath("C:\Users\u0117545\OneDrive - KU Leuven\2.Dataprocessing\Matlab\addons")

Timepoint   = 'T0';
movement    = "F";
path.root   = 'C:\Users\u0117545\KU Leuven\An De Groef - DATA';
path.out    = fullfile(path.root,'Output','Database_MovQual.mat');
plot_or_not = 1;

Affected_table = readtable(fullfile(path.root,"Aangedane zijde.xlsx"));



%% 2. load data
for subj = 1:12
    if subj < 10
        subj_name   = ['BC_00' num2str(subj)];
    elseif subj < 100
        subj_name   = ['BC_0' num2str(subj)];
    else
        subj_name   = ['BC_', num2str(subj)];
    end

    affected = Affected_table(strcmp(Affected_table.ppID, subj_name), "involved");


    disp(' ')
    disp(['Processing ' subj_name ': ' Timepoint '.....'])

    path.subj   = fullfile(path.root, subj_name, 'Xsens', Timepoint, 'Reproces');
    check_subj  = exist(path.subj);

    if check_subj == 7
        %initialize counters
        counterR        = 0;
        counterR_SSS    = 0;
        counterL        = 0;
        counterL_SSS    = 0;

        content = dir(path.subj);
        nfiles = size(content,1);

        % Start loop through ULIFT files per subject
        for file = 1:nfiles
            if contains(content(file).name, movement) && contains(content(file).name, '.mvnx') && ~contains(content(file).name, 'AF') && ~contains(content(file).name, 'ULIFT')
                file_ik = fullfile(path.subj, content(file).name);

                [~,name, ~] = fileparts(content(file).name);
                [fileName] = regexprep(name, '-', '_');

                d = strfind(name,'_');
                arm = content(file).name(d+1);

                if strcmp(arm, 'L')
                    sensorno    = 10;
                    segmentno   = 14;
                    jointno     = 12;

                    if strcmp(Affected_table{strcmp(Affected_table.ppID, subj_name), "involved"}, 'L')
                        side = 'affected';
                    else
                        side = 'unaffected';

                    end

                elseif strcmp(arm, 'R')
                    sensorno    = 6;
                    segmentno   = 10;
                    jointno     = 8;

                    if strcmp(Affected_table{strcmp(Affected_table.ppID, subj_name), "involved"}, 'R')
                        side = 'affected';
                    else
                        side = 'unaffected';
                    end
                end

                disp(['     ' 'Analysing: ' fileName '.....'])
                disp(['   ' 'Arm of interst: ' arm '.....'])

                %% 2.1 Load xsens data
                % Change the filename here to the name of the file you would like to import
                disp(['    ' content(file).name ': read xsens file'])
                [sensorData, segmentData, jointData]= MVN(file_ik);

                x = sensorData(sensorno).sensorFreeAcceleration(:,1);
                y = sensorData(sensorno).sensorFreeAcceleration(:,2);
                z = sensorData(sensorno).sensorFreeAcceleration(:,3);
                res = vecnorm(sensorData(sensorno).sensorFreeAcceleration,2,2);
                acc = table(x, y, z, res);
% 
%                 figure;
%                 stackedplot(acc)



                %% event detection || seperation of the different repetitions
                disp(['    ' content(file).name ': define seperate repetitions'])

                %filter data
                fc = 2;  %cutoff freq
                fs = 60; %sample freq
                [b,a] = butter(2, fc/(fs/2));

                velocity = filtfilt(b,a, segmentData(segmentno).velocity);
                velocityX = velocity(:,1);
                velocityY = velocity(:,2);
                velocityZ = velocity(:,3);
                velocityVec = vecnorm(velocity, 2,2);

                SensorFree = filtfilt(b,a, sensorData(sensorno).sensorFreeAcceleration);
                SensorFreeX = SensorFree(:,1);
                SensorFreeY = SensorFree(:,2);
                SensorFreeZ = SensorFree(:,3);
                SensorFreeVec = vecnorm(SensorFree,2,2);
                SensorFreeDiff = [diff(SensorFreeVec); 0];

                angularVel_LA = filtfilt(b,a, segmentData(segmentno).angularVelocity);
                angularVelX = angularVel_LA(:,1);
                angularVelY = angularVel_LA(:,2);
                angularVelZ = angularVel_LA(:,3);
                angularVelVec = vecnorm(angularVel_LA, 2, 2);
                angularVelDiff = [diff(angularVelVec); 0];

                % dataframes
                df.vel      = table(velocityX, velocityY, velocityZ, velocityVec);
                df.SenAcc   = table(SensorFreeX, SensorFreeY, SensorFreeZ, SensorFreeVec, SensorFreeDiff);
                df.Avel     = table(angularVelX, angularVelY, angularVelZ, angularVelVec, angularVelDiff);

                clear velocity velocityX velocityY velocityZ velocityVec
                clear sensorFree sensorFreeX sensorFreeY SensorFreeZ sensorFreeVec SensorFreeDiff
                clear angularVel_X angularVelY angularVelZ angularVelVec angularVelDiff


%                 figure; 
%                 stackedplot(df.Avel)
% 
%                 figure; 
%                 stackedplot(df.SenAcc)

                figure;
                stackedplot(df.vel)
                title([subj_name, ' ', fileName])
                %% LDLJ_A



                %% sample entropy
                %----------------
                r = 0.2;
                m = 2;
                SampleEntropy_x = sampen(acc.x, m, r);
                SampleEntropy_y = sampen(acc.y, m, r);
                SampleEntropy_z = sampen(acc.z, m, r);
                SampleEntropy_res = sampen(acc.res, m, r);


                ppID = string(subj_name);
                trial = string(fileName);

                if strcmp(side, 'affected')
                    SampEn_aff(subj, :) = table(ppID, trial, SampleEntropy_x, SampleEntropy_y, SampleEntropy_z, SampleEntropy_res);

                elseif strcmp(side, 'unaffected')
                    SampEn_unaff(subj, :) = table(ppID, trial, SampleEntropy_x, SampleEntropy_y, SampleEntropy_z, SampleEntropy_res);
                end

            end% file name contains movement and .mvnx
            %% structure the data for appending to the table

        end% loop though the number of files
    end% check if data folder exists
end% loop through number of subjects