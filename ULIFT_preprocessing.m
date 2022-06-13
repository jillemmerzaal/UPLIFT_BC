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
path.table  = fullfile(path.root,'Output');
plot_or_not = 1;

%% 2. load data
for subj = 3
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
        %T = readtable(fullfile(path.root, 'Output', [subj_name, '.xlsx']));

        %initialize counters
        counterR        = 0;
        counterR_SSS    = 0;
        counterL        = 0;
        counterL_SSS    = 0;
        
        counterRUN      = 0;
        content = dir(path.subj);
        nfiles = size(content,1);

        % Start loop through ULIFT files per subject
        for file = 1:nfiles
            if contains(content(file).name, movement) && contains(content(file).name, '.mvnx')
                number  = str2num(content(file).name(13:end-5));
                file_ik = fullfile(path.subj, content(file).name);

                [~,name, ~] = fileparts(content(file).name);
                [fileName] = regexprep(name, '-', '_');

%                 idx = find(ismember(T.filename, fileName));

                % if T.run(idx)==1

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

                %% 2.1 Load xsens data
                % Change the filename here to the name of the file you would like to import
                disp(['    ' content(file).name ': read xsens file'])
                [sensorData, segmentData, jointData]= MVN(file_ik);

                if contains(arm, 'L')
                    jointno     = 14;
                    segmentno   = 14;
                    sensorno    = 10;
                else
                    jointno     = 10;
                    segmentno   = 10;
                    sensorno    = 6;
                end


                %% 2.2 Define data needed for segmentation
                %-----------------------------------------
                disp(['    ' content(file).name ': define start and end points'])

                %filter data
                fc = 2;  %cutoff freq
                fs = 60; %sample freq
                [b,a] = butter(2, fc/(fs/2));

                position = filtfilt(b,a, segmentData(segmentno).position);
                positionX = position(:,1);
                positionY = position(:,2);
                positionZ = position(:,3);
                positionVec = vecnorm(position, 2,2);

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
                df.pos      = table(positionX, positionY, positionZ, positionVec);
                df.SenAcc   = table(SensorFreeX, SensorFreeY, SensorFreeZ, SensorFreeVec, SensorFreeDiff);
                df.Avel     = table(angularVelX, angularVelY, angularVelZ, angularVelVec, angularVelDiff);

                clear position positionX positionY positionZ positionVec
                clear sensorFree sensorFreeX sensorFreeY SensorFreeZ sensorFreeVec SensorFreeDiff
                clear angularVel_X angularVelY angularVelZ angularVelVec angularVelDiff

                %% set counters
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

                counterRUN = counterRUN + 1;
                %% Step 1: define change points based on position data
                %-----------------------------------------------------
                [changeIndices,segmentMean] = ischange(df.pos.positionZ,"MaxNumChanges",2);
                x = find(changeIndices);


                % check the number of highest points -- should be 3 per phase
                % If there are more than 3, trial will not be run.  
                %------------------------------------------------------------
                thresh_ph1 = mean(segmentMean(1:x(1))) + mean(segmentMean(1:x(1)))*0.05;
                [maxIndices_ph1, peakMag] = peakfinder(df.pos.positionZ(1:x(1)), [], thresh_ph1, 1, []);

                if length(maxIndices_ph1) > 3
                    thresh_ph1 = mean(peakMag) - 2*std(peakMag);
                    [maxIndices_ph1, ~] = peakfinder(df.pos.positionZ(1:x(1)), [], thresh_ph1, 1, []);
                end

                


                thresh_ph4 = mean(segmentMean(x(2):end)) + mean(segmentMean(x(2):end))*0.025;

                [maxIndices_ph4, peakMag] = peakfinder(df.pos.positionZ(x(2):end), [], thresh_ph4);

                if length(maxIndices_ph4) > 3
                    thresh_ph4 = mean(peakMag) - std(peakMag);
                    [maxIndices_ph4, ~] = peakfinder(df.pos.positionZ(x(2):end), [], thresh_ph4, 1, []);
                end

                maxIndices_ph4 = maxIndices_ph4 + x(2);

                ppID{counterRUN, 1} = subj_name;
                Phase1(counterRUN,1) = size(maxIndices_ph1,1);
                filename{counterRUN,1} = fileName;
                Phase4(counterRUN,1) = size(maxIndices_ph4,1);
                Tpoint{counterRUN, 1} = Timepoint;


                if Phase1(counterRUN,1) == 3 && Phase4(counterRUN,1) == 3
                    run(counterRUN,1) = 1;
                else
                    run(counterRUN,1 ) = 0;
                end
               

                if run(counterRUN,1) ==1
                    %% Step 2: Start and end of the ULIFT task
                    % define start and end of the ULIFT task as using the
                    % peak rate of thange of the angular velocity vector of
                    % the lower arm sensor
                    %------------------------------------------------------

                    % Start phase 1 
                    [peakLoc, peakMag] = peakfinder(df.Avel.angularVelDiff(1:x(1)));
                    localmax.all = peakLoc;
                    Thresh = mean(peakMag) *1.5;
                    [peakLoc, peakMag] = peakfinder(df.Avel.angularVelDiff(1:x(1)), [], Thresh);
                    localmax.thresh = peakLoc;

                    if isempty(localmax.thresh)
                        startPhase1 = localmax.all(1);
                    else
                        startPhase1 = localmax.thresh(1);
                    end

                    clear localmax peakLoc peakMag

                    % End phase 4
                    temp = df.Avel.angularVelDiff(x(2):end)*-1;
                    [peakLoc, peakMag] = peakfinder(temp);
                    localmin.all = peakLoc + x(2);

                    N = 1:height(df.Avel.angularVelDiff);
                    % localmin.all = min;
                    average = mean(peakMag);
                    Thresh = average *1.5;
                    [peakLoc, peakMag] = peakfinder(temp, [], Thresh);

                    localmin.thresh = peakLoc + x(2);
                    localmin.prominent = N(localmin.thresh)+x(2);
                    localmin.incon = N(localmin.all) + x(2);

                    if isempty(localmin.thresh)
                        Thresh = average ;
                        [peakLoc, peakMag] = peakfinder(temp, [], Thresh);
                        endPhase4 = peakLoc(end) + x(2);
                    else
                        endPhase4 = peakLoc(end) + x(2);
                    end

                    clear localmin peakLoc peakMag temp

                    %% Step 3: End phase 1 and start phase 4 
                    % define the end of phase 1 and the start of phase 4
                    % using the peak acceleration signal in x-direction of 
                    % thelower arm sensor
                    %------------------------------------------------------

                    % End phase 1
                    [temp, P] = islocalmin(df.SenAcc.SensorFreeX(1:x(1)));

                    localmin.all = temp;
                    Thresh = mean(P(localmin.all));

                    clear temp P
                    [temp, P] = islocalmin(df.SenAcc.SensorFreeX(1:x(1)), 'MinProminence',Thresh);
                    localmin.thresh = temp;
                    N = 1:height(df.SenAcc.SensorFreeX);

                    %select the less prominent minima between the last two most prominent minima
                    localmin.prominent = N(localmin.thresh);
                    localmin.incon = N(localmin.all);

                    % endPhase1_new = localmin.incon(find(localmin.incon == localmin.prominent(end))-1)
                    endPhase1 = localmin.incon(end-1);
                    clear localmin P temp

                    % start phase 4
                    temp = df.SenAcc.SensorFreeX(x(2):end);
                    [minima, P] = islocalmin(temp);

                    localmin.all = minima;
                    Thresh = mean(P(localmin.all)) + std(P(localmin.all)) *0.25;
                    clear minima P
                    [minima, P] = islocalmin(temp, 'MinProminence',Thresh);
                    localmin.thresh = minima;

                    % the first prominent acceleration peak
                    localmin.prominent = N(localmin.thresh)+ x(2);
                    localmin.incon = N(localmin.all)+ x(2);
                    startPhase4 = localmin.prominent(1);

                    clear localmin temp P
                   
                    % the time axis of the phases. 
                    T_phase1 = startPhase1:endPhase1;
                    T_phase4 = startPhase4:endPhase4;

                    if size(T_phase1,2) > 50 && size(T_phase4, 2) > 50
                        %% Step 4: Extract the relevant kinematics
                        %-------------------------------------
                        disp(['    ' content(file).name ': extract relevant Kinematics'])

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

                        for jnt = 1:4
                            % initiation of joint names
                            IK_X = [jointNames{jnt}, '_abbuction'];
                            IK_Y = [jointNames{jnt}, '_rotation'];
                            IK_Z = [jointNames{jnt}, '_flexion'];


                            %% Phase 1
                            %---------
                            temp.X_phase1 = jointData(jointNo(jnt)).jointAngle(T_phase1, 1); %fist phase
                            temp.Y_phase1 = jointData(jointNo(jnt)).jointAngle(T_phase1, 2); %fist phase
                            temp.Z_phase1 = jointData(jointNo(jnt)).jointAngle(T_phase1, 3);

                            if size(temp.X_phase1, 1) < 100
                                nf = 102;
                            else
                                nf = 101;
                            end
                            % Time normalised phases; phase 1
                            %--------------------------------
                            Data_out.(movement).(Timepoint).IK.(arm).Phase1.normalised.(IK_X)(:,counter) = interp1([1:size(temp.X_phase1,1)],...
                                temp.X_phase1', [1:(size(temp.X_phase1,1))/nf:size(temp.X_phase1,1)], 'spline');

                            Data_out.(movement).(Timepoint).IK.(arm).Phase1.normalised.(IK_Y)(:,counter) = interp1([1:size(temp.Y_phase1,1)],...
                                temp.Y_phase1', [1:(size(temp.Y_phase1,1))/nf:size(temp.Y_phase1,1)], 'spline');

                            Data_out.(movement).(Timepoint).IK.(arm).Phase1.normalised.(IK_Z)(:,counter) = interp1([1:size(temp.Z_phase1,1)],...
                                temp.Z_phase1', [1:(size(temp.Z_phase1,1))/nf:size(temp.Z_phase1,1)], 'spline');


                            %% phase 4
                            %---------

                            temp.X_phase4 = jointData(jointNo(jnt)).jointAngle(T_phase4, 1); %fourth phase
                            temp.Y_phase4 = jointData(jointNo(jnt)).jointAngle(T_phase4, 2); %fourth phase
                            temp.Z_phase4 = jointData(jointNo(jnt)).jointAngle(T_phase4, 3); %fourth phase

                            if size(temp.X_phase4, 1) < 100
                                nf = 102;
                            else
                                nf = 101;
                            end

                            % Time normalised phases; phase 4
                            %--------------------------------
                            Data_out.(movement).(Timepoint).IK.(arm).Phase4.normalised.(IK_X)(:,counter) = interp1([1:size(temp.X_phase4,1)],...
                                temp.X_phase4', [1:(size(temp.X_phase4,1))/nf:size(temp.X_phase4,1)], 'spline');

                            Data_out.(movement).(Timepoint).IK.(arm).Phase4.normalised.(IK_Y)(:,counter) = interp1([1:size(temp.Y_phase4,1)],...
                                temp.Y_phase4', [1:(size(temp.Y_phase4,1))/nf:size(temp.Y_phase4,1)], 'spline');

                            Data_out.(movement).(Timepoint).IK.(arm).Phase4.normalised.(IK_Z)(:,counter) = interp1([1:size(temp.Z_phase4,1)],...
                                temp.Z_phase4', [1:(size(temp.Z_phase4,1))/nf:size(temp.Z_phase4,1)], 'spline');


                            %% non-normalised data
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


                            clear temp
                        end

                        %% Display the results
                        %---------------------
                        if plot_or_not
                            figure;
                            tiledlayout('flow')

                            % display the results of the change points
                            nexttile
                            plot(df.pos.positionZ,"Color",[77 190 238]/255,"DisplayName","Input data")
                            hold on

                            % Plot segments between change points
                            plot(segmentMean,"Color",[64 64 64]/255,"DisplayName","Segment mean")

                            %Plot change points
                            x_rep = repelem(find(changeIndices),3);
                            y = repmat([ylim(gca) missing]',nnz(changeIndices),1);
                            plot(x_rep,y,"Color",[51 160 44]/255,"LineWidth",1,"DisplayName","Change points")
                            title("Number of change points: " + nnz(changeIndices))

                            hold off
                            %legend('Position',[0.85,0.25,0.15,0.2])
                            clear segmentMean x_rep y posdata

                            % Plot the start and end points!
                            nexttile
                            plot(df.SenAcc.SensorFreeX ,"Color",[77 190 238]/255,"DisplayName","sensorfree acceleration")
                            hold on
                            xline(startPhase1, "Color", '#A2142F', "DisplayName",'StrPh1')
                            xline(endPhase1, "Color", '#A2142F', "DisplayName",'EndPh1')

                            xline(startPhase4, "Color", '#EDB120',"LineWidth",1, "DisplayName",'StrPh4')
                            xline(endPhase4, "Color", '#EDB120', "LineWidth",1,"DisplayName",'EndPh4')
                            hold off
                            %legend('Position',[0.85,0.25,0.15,0.2])
                            title("Start/end points")
                            ylabel("Acceleration X")
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
                        end % too few data after segmentation
                    end%plot_or_not
                end% first check
            end % end if movement && .mvnx
        end %end number of files

        if strcmp(Timepoint, 'T0')
            T = table(ppID, filename, Tpoint, Phase1, Phase4, run);
            writetable(T, fullfile(path.table, [subj_name, '.xlsx']))
        else
            tablefile = fullfile(path.table, [subj_name, '.xlsx']);
            T_orig = readtable(tablefile);

            T_new = table(ppID, filename, Tpoint, Phase1, Phase4, run);
            T = [T_orig; T_new];            
            writetable(T, fullfile(path.table, [subj_name, '.xlsx']))


        end

    end% end check if subject path exists

    %% Save data
    %--------------

        if exist('Data_out','var')
            if exist(path.out,'file')
                load(path.out)
            end
    
            [Data.(subj_name).ULIFT.(Timepoint)] = Data_out.ULIFT.(Timepoint);
            save(path.out,'Data')
            clear Data Data_out
        end

    disp(['*********Finished ' subj_name '**********'])
    disp(' ')

end %end subjects
