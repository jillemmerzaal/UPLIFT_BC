%% maximal range of motion task
% This code extracts the maximal range of motion as calculated by the xsens
% system. Data collected with Xsens MVN 2021.2.
% Nescesary functions: 
%       1) MVN.m 
%       2) load_mvnx.m
% code written by 
%       dr. Jill Emmerzaal
%       KU Leuven, Tervuursevest 101, box 1501
%       Research Group for Rehabilitation in Internal Disorders

clearvars; %close all; clc
%% 1. input data
cd("C:\Users\u0117545\Documents\GitHub\ULIFT_BC")
addpath("C:\Users\u0117545\OneDrive - KU Leuven\2.Dataprocessing\Matlab\addons")

Timepoint   = 'T0';
activities  =  {'ABD'};%, 'AF', 'EXO'}; %Abductie'; %Abductie Anteflexie Exorotatie
path.root   = 'C:\Users\u0117545\KU Leuven\An De Groef - DATA';
path.out    = fullfile(path.root,'Output','Database_ROM.mat');
plot_or_not = 1;

%% 2. load data
for subj = 1%:8
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
    check_subj  = exist(path.subj, 'dir');

    if check_subj == 7
        %initialize counters
        counterR = 0;
        counterL = 0;

        content = dir(path.subj);
        nfiles = size(content,1);

        % Start loop through ROM files per subject
        for mv = 1:3
            movement = string(activities(mv));

            for file = 1:nfiles
                if contains(content(file).name, movement) && contains(content(file).name, '.mvnx')
                    file_ik = fullfile(path.subj, content(file).name);

                    [~,name, ~] = fileparts(content(file).name);
                    [fileName] = regexprep(name, '-', '_');

                    d = strfind(name,'_');
                    arm = content(file).name(d+1);

                    %% 2.1 Load xsens data
                    % Change the filename here to the name of the file you would like to import
                    disp(['    ' content(file).name ': read xsens file'])
                    [sensorData, segmentData, jointData]= MVN(file_ik);

                    if strcmp(arm, 'L')
                        sensorno    = 10;
                        segmentno   = 14;
                        jointno     = 12;
                    else
                        sensorno    = 6;
                        segmentno   = 10;
                        jointno     = 8;
                    end

                    %% extract the relevant peak kinematics from the three repetitions from each activity
                    switch movement
                        case 'ABD'
                            disp(['    ' movement ': Extract peak kinematics'])

                            df = jointData(jointno).jointAngleXZY;
                            figure;
                            plot(df)

                            [peakLoc, peakMag] = peakfinder(df(:,1));

                            if strcmp(arm, 'L')
                                maxABD_L = mean(peakMag);
                            else
                                maxABD_R = mean(peakMag);
                            end

                        case 'AF'
                            disp(['    ' movement ': Extract peak kinematics'])
                            df = jointData(jointno).jointAngle;

                            figure;
                            plot(df)

                            [peakLoc, peakMag] = peakfinder(df(:,3));

                            if strcmp(arm, 'L')
                                maxAF_L = mean(peakMag);
                            else
                                maxAF_R = mean(peakMag);
                            end

                        case 'EXO'
                            disp(['    ' movement ': Extract peak kinematics'])

                            df = jointData(jointno).jointAngleXZY;

                            figure;
                            plot(df)

                            [peakLoc, peakMag] = peakfinder(df(:,2),[],[],-1);

                            if strcmp(arm, 'L')
                                maxEXO_L = abs(mean(peakMag));
                            else
                                maxEXO_R = abs(mean(peakMag));
                            end
                    end% switch case
                end% end if movement && .mvnx
            end  % loop of number of files
        end%loop though activities
        %% structure the data for appending to the table
        subj_id = {subj_name};
        Data_out = {subj_id, maxABD_L, maxABD_R, maxAF_L, maxAF_R, maxEXO_L, maxEXO_R};
        close all
    end %check if subject exists
    %% Save data
    %--------------

    if exist('Data_out','var')
        if exist(path.out,'file')
            load(path.out)

            [Data] = [Data;Data_out];

        else 
            Data = table(subj_id, maxABD_L, maxABD_R, maxAF_L, maxAF_R, maxEXO_L, maxEXO_R);
        end


        save(path.out,'Data')
        clear Data Data_out
    end

    disp(['*********Finished ' subj_name '**********'])
    disp(' ')
    
end %end number of subjects
