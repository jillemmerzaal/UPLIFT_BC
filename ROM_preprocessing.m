%% maximal range of motion task
%
clear all; close all; clc
%% 1. input data
cd("C:\Users\u0117545\Documents\GitHub\ULIFT_BC")
addpath("C:\Users\u0117545\OneDrive - KU Leuven\2.Dataprocessing\Matlab\addons")


Timepoint   = 'T0';
movement    =  'F'; %Abductie'; %Abductie Anteflexie Exorotatie
path.root   = 'C:\Users\u0117545\KU Leuven\An De Groef - DATA';
path.out    = fullfile(path.root,'Output','Database_ULIFT.mat');

plot_or_not = 1;



%% 2. load data
for subj = 1%1:3
    if subj < 10
        subj_name   = ['BC_00' num2str(subj)];
    elseif subj < 100
        subj_name   = ['BC_0' num2str(subj)];
    else
        subj_name   = ['BC_', num2str(subj)];
    end

    disp(' ')
    disp(['Processing ' subj_name ': ' Timepoint '.....'])

    path.subj   = fullfile(path.root, subj_name, 'Xsens', Timepoint, 'Reproces');
    check_subj  = exist(path.subj);

    %% validation
    timepoints = readtable("C:\Users\u0117545\Documents\GitHub\ULIFT_BC\ValidationROM.xlsx");
    timepoints = timepoints(strcmp(timepoints.subj_id, subj_name),:);
    if check_subj == 7
        %initialize counters
        counterR = 0;
        counterL = 0;

        content = dir(path.subj);
        nfiles = size(content,1);

        % Start loop through ROM files per subject
        for file = 1:nfiles
            if contains(content(file).name, movement) && contains(content(file).name, '.mvnx')
                file_ik = fullfile(path.subj, content(file).name);

                [~,name, ~] = fileparts(content(file).name);
                [fileName] = regexprep(name, '-', '_');

                %                 if contains(content(file).name, 'R')
                %                     arm = 'right';
                %                 else
                %                     arm = 'left';
                %                 end



                d = strfind(name,'_');
                arm = content(file).name(d+1);

                interest = timepoints(strcmp(timepoints.Trial, name),:);


                

                %% 2.1 Load xsens data
                % Change the filename here to the name of the file you would like to import
                disp(['    ' content(file).name ': read xsens file'])
                [sensorData, segmentData, jointData]= MVN(file_ik);

                if contains(arm, 'L')
                    sensorno     = 10;
                    segmentno   = 14;
                else
                    sensorno     = 6;
                    segmentno   = 10;
                end
                %% Sensor and segment data
                %filter data
                fc = 2;  %cutoff freq
                fs = 60; %sample freq
                [b,a] = butter(2, fc/(fs/2));


                angularVel = filtfilt(b,a, segmentData(segmentno).angularVelocity);
                angularVelX = angularVel(:,1);
                angularVelY = angularVel(:,2);
                angularVelZ = angularVel(:,3);
                angularVelVec = vecnorm(angularVel, 2, 2);
                angularVelDiff = [diff(angularVelVec); 0];

                SensorFree = filtfilt(b,a, sensorData(sensorno).sensorFreeAcceleration);
                SensorFreeX = SensorFree(:,1);
                SensorFreeY = SensorFree(:,2);
                SensorFreeZ = SensorFree(:,3);
                SensorFreeVec = vecnorm(SensorFree,2,2);
                SensorFreeDiff = [diff(SensorFreeVec); 0];

                df.Avel     = table(angularVelX, angularVelY, angularVelZ, angularVelVec, angularVelDiff);
                df.SenAcc   = table(SensorFreeX, SensorFreeY, SensorFreeZ, SensorFreeVec, SensorFreeDiff);

                clear sensorFree sensorFreeX sensorFreeY SensorFreeZ sensorFreeVec SensorFreeDiff
                clear angularVel_X angularVelY angularVelZ angularVelVec angularVelDiff


                signals = fieldnames(df);
                for nPlot = 1:length(signals)
                    plottitle = {[signals{nPlot} ' data ' fileName]};

                    figure;
                    h = stackedplot(df.(signals{nPlot}));
                    title(plottitle)

                    %based on xsens
                    ax = findobj(h.NodeChildren, 'Type','Axes');
                    arrayfun(@(h)xline(h,interest.start_1,'LineWidth',1.5, "Color", '#A2142F', "DisplayName",'start1'),ax)
                    %arrayfun(@(h)xline(h,interest.start_2, 'LineWidth', 1.5,"Color", '#A2142F', "DisplayName",'start2'), ax)
                    arrayfun(@(h)xline(h,interest.start_3,'LineWidth',1.5, "Color", '#A2142F', "DisplayName",'start3'),ax)


                   
%                    arrayfun(@(h)xline(h,interest.stop_1,'LineWidth',1.5, "Color", 'blue', "DisplayName",'stop1'),ax)
%                     arrayfun(@(h)xline(h,interest.stop_2,'LineWidth',1.5, "Color", 'blue', "DisplayName",'stop1'),ax)
                    arrayfun(@(h)xline(h,interest.stop_3,'LineWidth',1.5, "Color", 'blue', "DisplayName",'stop1'),ax)

                end
            end % end if movement && .mvnx
        end % loop of number of files
    end %check if subject exists
    %% Save data
    %--------------

    if exist('Data_out','var')
        if exist(path.out,'file')
            load(path.out)
        end


        [Data.(subj_name).(movement).(Timepoint)] = Data_out.(movement).(Timepoint);
        save(path.out,'Data')
        clear Data Data_out
    end

    disp(['*********Finished ' subj_name '**********'])
    disp(' ')
end %end subjects
