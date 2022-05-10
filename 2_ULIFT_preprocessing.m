%% Extract the results of the ULIFT task
%

clear all; close all; clc
%% 1. input data
cd("C:\Users\u0117545\Documents\GitHub\ULIFT_BC")
addpath("C:\Users\u0117545\OneDrive - KU Leuven\2.Dataprocessing\Matlab\addons")


Timepoint   = 'T0';
movement    = "ULIFT";
path.root   = 'C:\Users\u0117545\KU Leuven\An De Groef - DATA';
path.out    = fullfile(path.root,'Output','Database_ULIFT.mat');


plot_or_not = 1;


%% 2. load data
for subj = 1
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
            if contains(content(file).name, movement) && contains(content(file).name, '.mvnx')
                number  = str2num(content(file).name(13:end-5));
                file_ik = fullfile(path.subj, content(file).name);

                [~,name, ~] = fileparts(content(file).name);
                [fileName] = regexprep(name, '-', '_');

                d = strfind(name,'_');
                if size(d,2) == 1
                    arm = content(file).name(d+1);

                elseif size(d,2) == 2
                    temp = content(file).name(d(2)+1);
                    arm = [temp, '_SSS']; % SSS = self-selected speed
                    clear temp
                end

                disp(['     ' 'Analysing: ' fileName '.....'])
                disp(['   ' 'Arm of interst: ' arm '.....'])


                %                 if contains(content(file).name, 'R')
                %                     arm = 'right';
                %                 else
                %                     arm = 'left';
                %                 end

                %% 2.1 Load xsens data
                % Change the filename here to the name of the file you would like to import
                disp(['    ' content(file).name ': read xsens file'])
                tree = load_mvnx(file_ik);

                % Read some basic data from the file
                mvnxVersion = tree.metaData.mvnx_version; % version of the MVN Studio used during recording

                if (isfield(tree.metaData, 'comment'))
                    fileComments = tree.metaData.comment; % comments written when saving the file
                end

                % Read some basic properties of the subject;

                frameRate = tree.metaData.subject_frameRate;
                suitLabel = tree.metaData.subject_label;
                originalFilename = tree.metaData.subject_originalFilename;
                recDate = tree.metaData.subject_recDate;
                segmentCount = tree.metaData.subject_segmentCount;

                % Retrieve sensor label
                % creates a struct with sensor data
                if isfield(tree,'sensorData') && isstruct(tree.sensorData)
                    sensorData = tree.sensorData;
                end

                % Retrieve segment labels
                % creates a struct with segment definitions
                if isfield(tree,'segmentData') && isstruct(tree.segmentData)
                    segmentData = tree.segmentData;
                end

                % Retrieve joint labels
                % creates a struct with segment definitions
                if isfield(tree,'jointData') && isstruct(tree.jointData)
                    jointData = tree.jointData;
                end

                if contains(arm, 'L')
                    jointno     = 14;
                    segmentno   = 14;
                else
                    jointno     = 10;
                    segmentno   = 10;
                end


                %% 2.2 Define start & end points of each repetition based on velocity of lower arm
                %---------------------------------------------------------------------------------
                disp(['    ' content(file).name ': define start and end points'])

                data = segmentData(jointno).velocity(:,3);

                %filter data
                fc = 2;  %cutoff freq
                fs = 60; %sample freq
                [b,a] = butter(2, fc/(fs/2));
                % freqz(b,a)
                datasmooth = filtfilt(b,a, data);

                %offset = mean(datasmooth)

                %datasmooth = datasmooth - mean(datasmooth);

                % change points
                % Find the change points of the position data and the minima and maxima in the
                % velocity data of the lower arm.
                %
                % Hypothesis is that on average the position of the high and the low sections
                % can be seperated using the position data. And a better idea of the precice start
                % end end points of the different target peaks
                %
                % Using both the changepoints and the arm velocity data a more robust segmenttation
                % of the ULIFT task can be achieved
                %----------------------------------

                peaks = segmentData(jointno).position(:,3);
                [changeIndices,segmentMean] = ischange(peaks,"MaxNumChanges",2);
                x = find(changeIndices);

                % Find local maxima and minima
                %-----------------------------
                n=1.5;
                thresh = mean(datasmooth) + n * std(datasmooth);
                [maxIndices, ~] = peakfinder(datasmooth, [], thresh, 1, []);
                [minIndices, ~] = peakfinder(datasmooth, [], thresh*-1, -1, []);


                % lower threshold untill 12 maxima are found. break while loop if the
                % threshold is lower than 0
                while size(maxIndices,1)<12
                    n=n-0.05;
                    thresh = mean(datasmooth) + n * std(datasmooth);
                    [maxIndices, ~] = peakfinder(datasmooth, [], thresh, 1, []);
                    [minIndices, ~] = peakfinder(datasmooth, [], thresh*-1, -1, []);
                    %disp(['threshold = ' num2str(thresh)])
                    if thresh < 0
                        break
                    end
                end

                % if while loop above had to break, and no 12 maxima are found,
                % try again but with a dicreasing data around the peaks
                idx = 4;
                sel = (max(datasmooth)-min(datasmooth))/idx;
                while size(maxIndices, 1) < 12
                    sel = sel - 0.05;
                    [maxIndices, ~] = peakfinder(datasmooth, sel, thresh, 1, []);
                    %disp(['Sel = ' num2str(sel)])

                end

                n=1.5;
                while size(minIndices,1) < 12
                    n=n-0.05;
                    thresh = (mean(datasmooth) + n * std(datasmooth));
                    [minIndices, ~] = peakfinder(datasmooth, [], thresh*-1, -1, []);


                    if thresh*-1 > -0.3
                        break
                    end

                end

                 idx = 4;
                sel = (max(datasmooth)-min(datasmooth))/idx;
                while size(minIndices, 1) < 12
                    sel = sel - 0.05;
                    [minIndices, ~] = peakfinder(datasmooth, sel, thresh*-1, -1, []);
                    %disp(['Sel = ' num2str(sel)])

                end

                % Determine start and end of phase 1 based on change points in the position
                % data, zero crossing and first 3 peaks in the velocity signal



                % start phase 1 --> highest shelf to middle shelf
                %------------zerocrossing vertical acceleration------------
                zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0.0001);                       % Returns Zero-Crossing Indices Of Argument Vector
                zx = zci(datasmooth);                                                           % Approximate Zero-Crossing Indices
                %find the crossings surounding the first positive peak
                %-----------------------------------------------------
                fistMax = maxIndices(1);
                startPhase1 = zx(find(zx < maxIndices(1),1,'last'));

                % end phase 1 --> highest shelf to middle shelf
                endPhase1_range = zx(find(zx < x(1), 100, "last"));

                % the end of phase one sould end right after the 3rd
                % minimum velocity peak, on the condition that we start
                % with a positive velocity peak.
                if maxIndices(1) < minIndices(1)
                    if sum(endPhase1_range > minIndices(3))
                        endPhase1 = endPhase1_range(find(endPhase1_range > minIndices(3), 1, 'first'));
                    else
                        endPhase1 = endPhase1_range(end);
                    end
                else
                    if sum(endPhase1_range > minIndices(4))
                        endPhase1 = endPhase1_range(find(endPhase1_range > minIndices(4), 1, 'first'));
                    else
                        endPhase1 = endPhase1_range(end);
                    end
                end


                % start and end of phase 4
                % Determine start and end points of phase 4 based on change points in the position
                % data, zero crossing and last 3 peaks in the velocity signal
                %------------------------------------------------------------
                endPhase4 = zx(find(zx > minIndices(end),1,'first'));

                startPhase4_range = zx(find(zx > x(2), 50, "first"));
                if sum(startPhase4_range > minIndices(end-3))
                    startPhase4 = startPhase4_range(find(startPhase4_range < minIndices(end-3), 1, 'last'));
                else
                    startPhase4 = startPhase4_range(1);
                end

                if isempty(startPhase4)
                    startPhase4 = zx(find(zx < maxIndices(end-3), 1, 'last'));
                end

                % start point of phase 4 should not be too far from the
                % third to last maximum peak. On the condition that we end
                % with a minimum peak.
                if maxIndices(end) > minIndices(end)
                    if abs(startPhase4 - maxIndices(end-3)) > 20
                        startPhase4 = startPhase4_range(find(startPhase4_range < maxIndices(end-3), 1, 'last'));
                    end
                else
                    if abs(startPhase4 - maxIndices(end-2)) > 20
                        startPhase4 = startPhase4_range(find(startPhase4_range < maxIndices(end-2), 1, 'last'));
                    end
                end


                adapted_range = startPhase4_range(find(startPhase4_range > maxIndices(end-3) ...
                    & startPhase4_range < maxIndices(end-2)));


                if size(adapted_range, 1) < 3
                    temp = adapted_range(1);
                else
                    temp = adapted_range(3);
                end

                
                if temp ~= startPhase4
                    startPhase4 = temp;
                end
                clear temp



                if isempty(startPhase1 )
                    startPhase1 = 1;
                end
                T_phase1 = startPhase1:endPhase1;
                T_phase4 = startPhase4:endPhase4;
                %% 2.3 Extract the relevant kinematics
                %-------------------------------------
                disp(['    ' content(file).name ': extract relevant Kinematics'])

                if strcmp(arm, 'R_SSS')
                    counterR_SSS    = counterR_SSS + 1;
                    counter         = counterR_SSS;
                    jointNo         = 7:10; % right upper extremity
                elseif strcmp(arm, 'R')
                    counterR        = counterR + 1;
                    counter         = counterR;
                    jointNo         = 7:10; % right upper extremity
                elseif strcmp(arm, 'L_SSS')
                    counterL_SSS    = counterL_SSS + 1;
                    counter         = counterL_SSS;
                    jointNo         = 11:14; % left upper extremity
                elseif strcmp(arm, 'L')
                    counterL        = counterL + 1;
                    counter         = counterL;
                    jointNo         = 11:14; % left upper extremity
                end

                % set General information per participant, per trial.
                Data_out.(movement).(Timepoint).General.CutIndices.(fileName) = [startPhase1, endPhase1, startPhase4, endPhase4];

                % initialise joint names
                jointNames = ['Scapula', "Glenohumeraal", "Elbow", "Wrist"];

                % Save the timecurves to struct, per participant, per
                % repetition.
                % Struct holds:
                % full timecurve
                % start and end points of phase 1 and phase 4
                % timecurves for each repetition
                % time normalised timecurves for each repetition
                %-------------------------------------------------
                nf = 101;
                for jnt = 1:4
                    % temporary time curves with the phase data
                    %------------------------------------------
                    temp.X_phase1 = jointData(jointNo(jnt)).jointAngle(T_phase1, 1); %fist phase
                    temp.X_phase4 = jointData(jointNo(jnt)).jointAngle(T_phase4, 1); %fourth phase
                    IK_X = [jointNames{jnt}, '_abbuction']; % initiation of joint names

                    temp.Y_phase1 = jointData(jointNo(jnt)).jointAngle(T_phase1, 2); %fist phase
                    temp.Y_phase4 = jointData(jointNo(jnt)).jointAngle(T_phase4, 2); %fourth phase
                    IK_Y = [jointNames{jnt}, '_rotation'];

                    temp.Z_phase1 = jointData(jointNo(jnt)).jointAngle(T_phase1, 3);
                    temp.Z_phase4 = jointData(jointNo(jnt)).jointAngle(T_phase4, 3); %fourth phase
                    IK_Z = [jointNames{jnt}, '_flexion'];

                    % Full ULIFT time data
                    %---------------------
                    Data_out.(movement).(Timepoint).IK.(arm).raw.(fileName).(IK_X) = jointData(jointNo(jnt)).jointAngle(:, 1);
                    Data_out.(movement).(Timepoint).IK.(arm).raw.(fileName).(IK_Y) = jointData(jointNo(jnt)).jointAngle(:, 2);
                    Data_out.(movement).(Timepoint).IK.(arm).raw.(fileName).(IK_Z) = jointData(jointNo(jnt)).jointAngle(:, 3);

                    % Timedata of the phases; phase 1
                    %---------------------------------
                    Data_out.(movement).(Timepoint).IK.(arm).Phase1.(fileName).(IK_X) = temp.X_phase1;
                    Data_out.(movement).(Timepoint).IK.(arm).Phase1.(fileName).(IK_Y) = temp.Y_phase1;
                    Data_out.(movement).(Timepoint).IK.(arm).Phase1.(fileName).(IK_Z) = temp.Z_phase1;

                    % Timedata of the phases; phase 4
                    %--------------------------------
                    Data_out.(movement).(Timepoint).IK.(arm).Phase4.(fileName).(IK_X) = temp.X_phase4;
                    Data_out.(movement).(Timepoint).IK.(arm).Phase4.(fileName).(IK_Y) = temp.Y_phase4;
                    Data_out.(movement).(Timepoint).IK.(arm).Phase4.(fileName).(IK_Z) = temp.Z_phase4;

                    % Time normalised phases; phase 1
                    %--------------------------------
                    Data_out.(movement).(Timepoint).IK.(arm).Phase1.normalised.(IK_X)(:,counter) = interp1([1:size(temp.X_phase1,1)],...
                        temp.X_phase1', [1:(size(temp.X_phase1,1))/nf:size(temp.X_phase1,1)], 'spline');

                    Data_out.(movement).(Timepoint).IK.(arm).Phase1.normalised.(IK_Y)(:,counter) = interp1([1:size(temp.Y_phase1,1)],...
                        temp.Y_phase1', [1:(size(temp.Y_phase1,1))/nf:size(temp.Y_phase1,1)], 'spline');

                    Data_out.(movement).(Timepoint).IK.(arm).Phase1.normalised.(IK_Z)(:,counter) = interp1([1:size(temp.Z_phase1,1)],...
                        temp.Z_phase1', [1:(size(temp.Z_phase1,1))/nf:size(temp.Z_phase1,1)], 'spline');

                    % Time normalised phases; phase 4
                    %--------------------------------
                    Data_out.(movement).(Timepoint).IK.(arm).Phase4.normalised.(IK_X)(:,counter) = interp1([1:size(temp.X_phase4,1)],...
                        temp.X_phase4', [1:(size(temp.X_phase4,1))/nf:size(temp.X_phase4,1)], 'spline');

                    Data_out.(movement).(Timepoint).IK.(arm).Phase4.normalised.(IK_Y)(:,counter) = interp1([1:size(temp.Y_phase4,1)],...
                        temp.Y_phase4', [1:(size(temp.Y_phase4,1))/nf:size(temp.Y_phase4,1)], 'spline');

                    Data_out.(movement).(Timepoint).IK.(arm).Phase4.normalised.(IK_Z)(:,counter) = interp1([1:size(temp.Z_phase4,1)],...
                        temp.Z_phase4', [1:(size(temp.Z_phase4,1))/nf:size(temp.Z_phase4,1)], 'spline');
                    clear temp
                end

                %% Display the results
                %---------------------
                if plot_or_not
                    figure;
                    tiledlayout('flow')

                    %display the results of the filtered and unfiltered data
                    nexttile
                    plot(data, "Color",[77 190 238]/255, 'DisplayName', "Unfiltered");
                    hold on;
                    plot(datasmooth, "Color",'#A2142F', "DisplayName","Filtered 1Hz")
                    hold off
                    title("filterd velocity data")

                    % display the results of the change points
                    nexttile
                    plot(peaks,"Color",[77 190 238]/255,"DisplayName","Input data")
                    hold on

                    % Plot segments between change points
                    plot(segmentMean,"Color",[64 64 64]/255,"DisplayName","Segment mean")

                    %Plot change points
                    x = repelem(find(changeIndices),3);
                    y = repmat([ylim(gca) missing]',nnz(changeIndices),1);
                    plot(x,y,"Color",[51 160 44]/255,"LineWidth",1,"DisplayName","Change points")
                    title("Number of change points: " + nnz(changeIndices))

                    hold off
                    %legend('Position',[0.85,0.25,0.15,0.2])
                    clear segmentMean x y peaks
                    %
                    % display the results of the peak detection
                    nexttile
                    plot(datasmooth,"Color",[77 190 238]/255,"DisplayName","Input data")
                    hold on
                    % Plot local maxima
                    plot(maxIndices,datasmooth(maxIndices),"^","Color",[217 83 25]/255,...
                        "MarkerFaceColor",[217 83 25]/255,"DisplayName","Local maxima")
                    % Plot local minima
                    plot(minIndices,datasmooth(minIndices),"v","Color",[237 177 32]/255,...
                        "MarkerFaceColor",[237 177 32]/255,"DisplayName","Local minima")
                    title("Number of extrema: " + (nnz(maxIndices)+nnz(minIndices)))
                    %legend('Position',[0.85,0.25,0.15,0.2])
                    % input the change points
                    x = find(changeIndices);
                    xline(x(1), "Color",[51 160 44]/255,"LineWidth",1, "DisplayName", "+/- endPh1")
                    xline(x(2), "Color",[51 160 44]/255,"LineWidth",1, "DisplayName", "+/- strPh4")

                    yline(0, "Color",[51 160 44]/255,"LineWidth",1, "DisplayName", "+/- zerocros")
                    hold off
                    %
                    % Plot the start and end points!
                    nexttile
                    plot(datasmooth,"Color",[77 190 238]/255,"DisplayName","Velocity")
                    hold on
                    xline(startPhase1, "Color", '#A2142F', "DisplayName",'StrPh1')
                    xline(endPhase1, "Color", '#A2142F', "DisplayName",'EndPh1')

                    xline(startPhase4, "Color", '#EDB120',"LineWidth",1, "DisplayName",'StrPh4')
                    xline(endPhase4, "Color", '#EDB120', "LineWidth",1,"DisplayName",'EndPh4')
                    hold off
                    %legend('Position',[0.85,0.25,0.15,0.2])
                    title("Start/end points")
                    ylabel("Velocity Z")
                    %
                    nexttile
                    plot(segmentData(jointno).position(:,3),"Color",[77 190 238]/255, "DisplayName", "position")
                    hold on
                    xline(startPhase1, "Color", '#A2142F', "DisplayName",'StrPh1')
                    xline(endPhase1, "Color", '#A2142F', "DisplayName",'EndPh1')

                    xline(startPhase4, "Color", '#EDB120',"LineWidth",1, "DisplayName",'StrPh4')
                    xline(endPhase4, "Color", '#EDB120', "LineWidth",1,"DisplayName",'EndPh4')
                    hold off
                    %legend('Position',[0.85,0.25,0.15,0.2])
                    title("Start/end points")
                    ylabel("position Z")

                    nexttile
                    stackedplot(Data_out.(movement).(Timepoint).IK.(arm).Phase1.normalised.Glenohumeraal_flexion)
                    title('Flexion/extension Shoulder--phase 1')

                    nexttile
                    stackedplot(Data_out.(movement).(Timepoint).IK.(arm).Phase4.normalised.Glenohumeraal_flexion);
                    title('Flexion/extension Shoulder--phase 4')

                    disp('      ')
                end
            end % end if movement && .mvnx
        end %end number of files
    end% end check if subject path exists

    %% Save data
    %--------------

    if exist('Data_out','var')
        if exist(path.out,'file')
            load(path.out)
        end


        subj_name2 = [subj_name, '_Reprocess'];
        [Data.(subj_name).ULIFT.(Timepoint)] = Data_out.ULIFT.(Timepoint);
        save(path.out,'Data')
        clear Data Data_out
    end

    disp(['*********Finished ' subj_name '**********'])
    disp(' ')

end %end subjects
