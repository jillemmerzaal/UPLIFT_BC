%% Extract the results of the ULIFT task
%

clear all; close all; clc
cd("C:\Users\u0117545\Documents\GitHub\ULIFT_BC")
addpath("C:\Users\u0117545\OneDrive - KU Leuven\2.Dataprocessing\Matlab\addons")


Timepoint   = 'T0';
movement    = "ULIFT";
path.root   = 'C:\Users\u0117545\KU Leuven\An De Groef - UPLIFT-BC\INVESTIGATOR SITE FILE\5. Data';
plot_or_not = 1;


%% Manual input is no longer nescesary after this point

% Start loop through subjects
for subj = 1
    if subj < 10
        subj_name   = ['BCT_00' num2str(subj)];
    elseif subj < 100
        subj_name   = ['BCT_0' num2str(subj)];
    else
        subj_name   = ['BCT_', num2str(subj)];
    end

    disp(' ')
    disp(['Processing ' subj_name ': ' Timepoint '.....'])


    path.subj   = fullfile(path.root, subj_name, 'Xsens', Timepoint);
    check_subj  = exist(path.subj);

    if check_subj == 7
        counterR = 0;
        counterL = 0;

        content = dir(path.subj);
        nfiles = size(content,1);

        % Start loop through ULIFT files per subject

        for file = 4%1:nfiles
            if contains(content(file).name, movement) && contains(content(file).name, '.mvnx')
                number  = str2num(content(file).name(13:end-5));
                file_ik = fullfile(path.subj, content(file).name);




                [~,name, ~] = fileparts(content(file).name);
                [fileName] = regexprep(name, '-', '_');


                if contains(content(file).name, 'R')
                    arm = 'right';
                else
                    arm = 'left';
                end


                % Change the filename here to the name of the file you would like to import
                tree = load_mvnx(file_ik);

                %% Read some basic data from the file
                mvnxVersion = tree.metaData.mvnx_version; % version of the MVN Studio used during recording

                if (isfield(tree.metaData, 'comment'))
                    fileComments = tree.metaData.comment; % comments written when saving the file
                end
                %%
                % Read some basic properties of the subject;

                frameRate = tree.metaData.subject_frameRate;
                suitLabel = tree.metaData.subject_label;
                originalFilename = tree.metaData.subject_originalFilename;
                recDate = tree.metaData.subject_recDate;
                segmentCount = tree.metaData.subject_segmentCount;


                %%
                % Retrieve sensor label

                %creates a struct with sensor data
                if isfield(tree,'sensorData') && isstruct(tree.sensorData)
                    sensorData = tree.sensorData;
                end
                %%
                % Retrieve segment labels

                %creates a struct with segment definitions
                if isfield(tree,'segmentData') && isstruct(tree.segmentData)
                    segmentData = tree.segmentData;
                end
                %%
                % Retrieve joint labels

                %creates a struct with segment definitions
                if isfield(tree,'jointData') && isstruct(tree.jointData)
                    jointData = tree.jointData;
                end

                if strcmp(arm, 'left')
                    jointno     = 14;
                    segmentno   = 14;
                else
                    jointno     = 10;
                    segmentno   = 10;
                end


                %% Define start & end points of each repetition based on velocity of lower arm

                data = segmentData(jointno).velocity(:,3);
                %filter data
                fc = 1;  %cutoff freq
                fs = 60; %sample freq
                [b,a] = butter(2, fc/(fs/2));
                % freqz(b,a)
                datasmooth = filtfilt(b,a, data);


                %% change points
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

                while size(maxIndices,1)<12
                    n=n-0.05;
                    thresh = mean(datasmooth) + n * std(datasmooth);
                    [maxIndices, ~] = peakfinder(datasmooth, [], thresh, 1, []);
                    [minIndices, ~] = peakfinder(datasmooth, [], thresh*-1, -1, []);

                end

                while size(minIndices,1) < 12
                    n=n-0.05;
                    thresh = mean(datasmooth) + n * std(datasmooth);
                    [minIndices, ~] = peakfinder(datasmooth, [], thresh*-1, -1, []);
                end

                %% Start and end of phase 1
                % Determine start and end of phase 1 based on change points in the position
                % data, zero crossing and first 3 peaks in the velocity signal

                % start phase 1 --> highest shelf to middle shelf
                %------------zerocrossing vertical acceleration------------
                zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0.0002);                       % Returns Zero-Crossing Indices Of Argument Vector
                zx = zci(datasmooth);                                                           % Approximate Zero-Crossing Indices

                %find the crossings surounding the first positive peak
                %-----------------------------------------------------
                fistMax = maxIndices(1);
                startPhase1 = zx(find(zx < maxIndices(1),1,'last'));

                % end phase 1 --> highest shelf to middle shelf
                endPhase1 = zx(find(zx < x(1), 30, "last"));

                if sum(endPhase1 > maxIndices(4))
                    endPhase1 = endPhase1(find(endPhase1 < maxIndices(4), 1, 'last'));
                else
                    endPhase1 = endPhase1(end);
                end

                %% start and end of phase 4
                % Determine start and end points of phase 4 based on change points in the position
                % data, zero crossing and last 3 peaks in the velocity signal
                %------------------------------------------------------------
                endPhase4 = zx(find(zx > minIndices(end),1,'first'));

                startPhase4 = zx(find(zx > x(2), 10, "first"));
                if sum(startPhase4 > minIndices(end-3))
                    startPhase4 = startPhase4(find(startPhase4 < minIndices(end-3), 1, 'last'));
                else
                    startPhase4 = startPhase4(1);
                end

                if isempty(startPhase4)
                    startPhase4 = zx(find(zx < minIndices(end-3), 1, 'last'));
                end

                %% Extract the relevant timeseries
                %---------------------------------

                if strcmp(arm, 'right')
                    counterR    = counterR + 1;
                    counter     = counterR;
                else
                    counterL    = counterL + 1;
                    counter     = counterL;
                end


                T_phase1 = startPhase1:endPhase1;
                T_phase4 = startPhase4:endPhase4;

                if strcmp(arm, 'left')
                    jointNo     = 11:14; % left upper extremity
                    segmentNo   = 12:15;
                else
                    jointNo     = 7:10; % right upper extremity
                    segmentNo   = 8:11;
                end

                jointNames = ['Scapula', "Glenohumeraal", "Elbow", "Wrist"];

                % timenormalise everything to 100% aka 101 samples
                %-------------------------------------------------
                nf = 101;
                for jnt = 1:4
                    temp.X_phase1 = jointData(jointNo(jnt)).jointAngle(T_phase1, 1); %fist phase
                    temp.X_phase4 = jointData(jointNo(jnt)).jointAngle(T_phase4, 1); %fourth phase
                    IK_X = [jointNames{jnt}, 'X'];

                    temp.Y_phase1 = jointData(jointNo(jnt)).jointAngle(T_phase1, 2);
                    temp.Y_phase4 = jointData(jointNo(jnt)).jointAngle(T_phase4, 2); %fourth phase

                    IK_Y = [jointNames{jnt}, 'Y'];

                    temp.Z_phase1 = jointData(jointNo(jnt)).jointAngle(T_phase1, 3);
                    temp.Z_phase4 = jointData(jointNo(jnt)).jointAngle(T_phase4, 3); %fourth phase

                    IK_Z = [jointNames{jnt}, 'Z'];

                    %normalize data to 101 frames
                    %----------------------------
                    if strcmp(arm, 'left')
                        subj_id = [subj_name, 'L'];

                        % Timedata of the individual phases; phase 1
                        %-------------------------------------------
                        Kinematics.(subj_id).(IK_X).Phase1.(fileName) = temp.X_phase1;
                        Kinematics.(subj_id).(IK_Y).Phase1.(fileName) = temp.Y_phase1;
                        Kinematics.(subj_id).(IK_Z).Phase1.(fileName) = temp.Z_phase1;

                        % Timedata of the individual phases; phase 4
                        %-------------------------------------------
                        Kinematics.(subj_id).(IK_X).Phase4.(fileName) = temp.X_phase4;
                        Kinematics.(subj_id).(IK_Y).Phase4.(fileName) = temp.Y_phase4;
                        Kinematics.(subj_id).(IK_Z).Phase4.(fileName) = temp.Z_phase4;

                        % Time normalised phases; phase 1
                        %--------------------------------
                        Kinematics.(subj_id).(IK_X).Phase1.normalised(:,counter) = interp1([1:size(temp.X_phase1,1)],...
                            temp.X_phase1', [1:(size(temp.X_phase1,1))/nf:size(temp.X_phase1,1)], 'spline');

                        Kinematics.(subj_id).(IK_Y).Phase1.normalised(:,counter) = interp1([1:size(temp.Y_phase1,1)],...
                            temp.Y_phase1', [1:(size(temp.Y_phase1,1))/nf:size(temp.Y_phase1,1)], 'spline');

                        Kinematics.(subj_id).(IK_Z).Phase1.normalised(:,counter) = interp1([1:size(temp.Z_phase1,1)],...
                            temp.Z_phase1', [1:(size(temp.Z_phase1,1))/nf:size(temp.Z_phase1,1)], 'spline');

                        % Time normalised phases; phase 4
                        %--------------------------------
                        Kinematics.(subj_id).(IK_X).Phase4.normalised(:,counter) = interp1([1:size(temp.X_phase4,1)],...
                            temp.X_phase4', [1:(size(temp.X_phase4,1))/nf:size(temp.X_phase4,1)], 'spline');

                        Kinematics.(subj_id).(IK_Y).Phase4.normalised(:,counter) = interp1([1:size(temp.Y_phase4,1)],...
                            temp.Y_phase4', [1:(size(temp.Y_phase4,1))/nf:size(temp.Y_phase4,1)], 'spline');

                        Kinematics.(subj_id).(IK_Z).Phase4.normalised(:,counter) = interp1([1:size(temp.Z_phase4,1)],...
                            temp.Z_phase4', [1:(size(temp.Z_phase4,1))/nf:size(temp.Z_phase4,1)], 'spline');


                    else
                        subj_id = [subj_name, 'R'];

                        % Timedata of the individual phases; phase 1
                        %-------------------------------------------
                        Kinematics.(subj_id).(IK_X).Phase1.(fileName) = temp.X_phase1;
                        Kinematics.(subj_id).(IK_Y).Phase1.(fileName) = temp.Y_phase1;
                        Kinematics.(subj_id).(IK_Z).Phase1.(fileName) = temp.Z_phase1;

                        % Timedata of the individual phases; phase 4
                        %-------------------------------------------
                        Kinematics.(subj_id).(IK_X).Phase4.(fileName) = temp.X_phase4;
                        Kinematics.(subj_id).(IK_Y).Phase4.(fileName) = temp.Y_phase4;
                        Kinematics.(subj_id).(IK_Z).Phase4.(fileName) = temp.Z_phase4;

                        % Time normalised phases; phase 1
                        %--------------------------------
                        Kinematics.(subj_id).(IK_X).Phase1.normalised(:,counter) = interp1([1:size(temp.X_phase1,1)],...
                            temp.X_phase1', [1:(size(temp.X_phase1,1))/nf:size(temp.X_phase1,1)], 'spline');

                        Kinematics.(subj_id).(IK_Y).Phase1.normalised(:,counter) = interp1([1:size(temp.Y_phase1,1)],...
                            temp.Y_phase1', [1:(size(temp.Y_phase1,1))/nf:size(temp.Y_phase1,1)], 'spline');

                        Kinematics.(subj_id).(IK_Z).Phase1.normalised(:,counter) = interp1([1:size(temp.Z_phase1,1)],...
                            temp.Z_phase1', [1:(size(temp.Z_phase1,1))/nf:size(temp.Z_phase1,1)], 'spline');

                        % Time normalised phases; phase 4
                        %--------------------------------
                        Kinematics.(subj_id).(IK_X).Phase4.normalised(:,counter) = interp1([1:size(temp.X_phase4,1)],...
                            temp.X_phase4', [1:(size(temp.X_phase4,1))/nf:size(temp.X_phase4,1)], 'spline');

                        Kinematics.(subj_id).(IK_Y).Phase4.normalised(:,counter) = interp1([1:size(temp.Y_phase4,1)],...
                            temp.Y_phase4', [1:(size(temp.Y_phase4,1))/nf:size(temp.Y_phase4,1)], 'spline');

                        Kinematics.(subj_id).(IK_Z).Phase4.normalised(:,counter) = interp1([1:size(temp.Z_phase4,1)],...
                            temp.Z_phase4', [1:(size(temp.Z_phase4,1))/nf:size(temp.Z_phase4,1)], 'spline');
                    end
                end
                clear temp

                %% Display the results
                %---------------------
                if plot_or_not
                    figure;

                    % Plot position of lower arm segement
                    tiledlayout('flow')
                    nexttile;
                    plot(tree.segmentData(jointno).position)
                    xlabel('frames')
                    ylabel('Position in the global frame')
                    legend('x','y','z')
                    title ('Position of lower arm segment')

                    % Plot 3D displacement of segment 14
                    %figure('name','Position of first segment in 3D')
                    nexttile;
                    plot3(tree.segmentData(jointno).position(:,1),tree.segmentData(jointno).position(:,2),tree.segmentData(jointno).position(:,3));
                    xlabel('x')
                    ylabel('y')
                    zlabel('z')
                    title ('Displacement of lower arm in space')


                    % display the results of the filtered and unfiltered data
                    nexttile; plot(data, "Color",[77 190 238]/255, 'DisplayName', "Unfiltered");
                    hold on;
                    plot(datasmooth, "Color",'#A2142F', "DisplayName","Filtered 1Hz")
                    hold off
                    title(content(file).name)

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
                end




            end % end if movement && .mvnx
        end %end number of files
    end% end check if subject path exists
end %end subjects
