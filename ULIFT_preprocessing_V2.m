%% Extract the results of the ULIFT task
%

clear all; close all; clc
%% 1. input data
cd("C:\Users\u0117545\Documents\GitHub\ULIFT_BC")
addpath("C:\Users\u0117545\OneDrive - KU Leuven\2.Dataprocessing\Matlab\addons")
addpath("C:\Users\u0117545\Documents\GitHub\ULIFT_BC")
Timepoint   = 'T0';
movement    = "ULIFT";
path.root   = 'C:\Users\u0117545\KU Leuven\An De Groef - DATA';
path.out    = fullfile(path.root,'Output','Database_ULIFT.mat');
path.table  = fullfile(path.root,'Output');
plot_or_not = 1;


%% 2. load data
for subj = 7%1:20
    if subj < 10
        subj_name   = ['BC_00' num2str(subj)];
    elseif subj < 100
        subj_name   = ['BC_0' num2str(subj)];
    else
        subj_name   = ['BC_', num2str(subj)];
    end


    diary myDiaryFile
    disp(' ')
    cprintf('*text', 'Processing: %s at Timepoint: %s....... \n', subj_name, Timepoint)

    path.subj   = fullfile(path.root, subj_name, 'Xsens', Timepoint, 'Reproces');
    check_subj  = exist(path.subj);

    if check_subj == 7
        %T = readtable(fullfile(path.root, 'Output', [subj_name, '.xlsx']));

        cd(path.subj)

        fprintf('\t current directory changed: %s \n', path.subj)
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
            try
                % try  to run the following code, if an error occurs
                % nothing will happen and no data will be exported.
                % A file will be created with the subject name and
                % timepoint that caused an error.
                if contains(content(file).name, movement) && contains(content(file).name, '.mvnx')
                    number  = str2num(content(file).name(13:end-5));

                    file_ik = content(file).name;

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

                    cprintf('*text', '\t Analysing: %s \n', fileName)
                    fprintf('\t\t Arm of interest: %s \n',  arm)

                    %% 2.1 Load xsens data
                    % Change the filename here to the name of the file you would like to import
                    fprintf('\t \t %s: read xsens file \n', content(file).name)

                    [sensorData, segmentData, jointData, tree]= MVN(file_ik);



                    if isfield(tree, 'fileComments')
                        fprintf('\t \t \t Comment in file: %s \n', tree.fileComments)
                    end

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
                    fprintf('\t \t %s: define start and end points \n', content(file).name)

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


                    if run(counterRUN,1) == 1
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
                        %% Step 2: Start and end of the ULIFT task
                        % define start and end of the ULIFT task as using the
                        % peak rate of thange of the angular velocity vector of
                        % the lower arm sensor
                        %------------------------------------------------------

                        % Start phase 1
                        [peakLoc, peakMag] = peakfinder(df.Avel.angularVelDiff(1:x(1)));
                        localmax.all = peakLoc;
                        Thresh = mean(peakMag) * 1.5;
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
                        %localmin.incon = N(localmin.all) + x(2);

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

                        %% middle phase 1
                        fprintf('\t \t %s: Segment middle arm movement phase 1 \n', content(file).name)

                        [peakLoc.phase1, peakMag.phase1] = peakfinder(df.SenAcc.SensorFreeX(T_phase1), [],[],-1);

                        peakLoc.phase1 = peakLoc.phase1 + startPhase1;
                        grab_idx.phase1 = peakLoc.phase1(1:2:end);
                        release_idx.phase1 = peakLoc.phase1(2:2:end);

                        T_middle_phase1 = grab_idx.phase1(2):grab_idx.phase1(3);

                        % throw error if the number of acceleration peaks exceed 6
                        %=========================================================
                        if size(peakLoc.phase1, 1) >7
                            fprintf('\t \t %s: Please select peaks phase 1 \n', content(file).name)



                            figure('Units','normalized','Position',[0 0 1 1]);
                            subplot(2,1,1);
                            % display the results of the change points

                            plot(df.pos.positionZ,"Color",[77 190 238]/255,"DisplayName","Input data")
                            hold on

                            % Plot segments between change points
                            plot(segmentMean,"Color",[64 64 64]/255,"DisplayName","Segment mean")

                            %Plot change points
                            x_rep = repelem(find(changeIndices),3);
                            y = repmat([ylim(gca) missing]',nnz(changeIndices),1);
                            plot(x_rep,y,"Color",[51 160 44]/255,"LineWidth",1,"DisplayName","Change points")
                            title("Number of change points: " + nnz(changeIndices))

                            xline(startPhase1, "Color", '#A2142F', "DisplayName",'StrPh1')
                            xline(endPhase1, "Color", '#A2142F', "DisplayName",'EndPh1')

                            xline(startPhase4, "Color", '#EDB120',"LineWidth",1, "DisplayName",'StrPh4')
                            xline(endPhase4, "Color", '#EDB120', "LineWidth",1,"DisplayName",'EndPh4')
                            hold off

                            %phase 1
                            subplot(2,1,2)
                            yline_phase1 = ones(size(T_phase1)) * (min(peakMag.phase1) * 1.2);

                            plot(df.SenAcc.SensorFreeX); hold on
                            plot(T_phase1, yline_phase1, 'LineWidth',10, 'Color', '#097392')
                            plot(grab_idx.phase1, peakMag.phase1(1:2:end),'o', 'MarkerSize', 7.5, 'MarkerFaceColor', 	'#A2142F', 'MarkerEdgeColor', 	'#A2142F')
                            plot(release_idx.phase1, peakMag.phase1(2:2:end), 'o', 'MarkerSize', 7.5, 'MarkerEdgeColor', 	'#0072BD', 'MarkerFaceColor', 	'#0072BD')

                            yline_middle_phase1 = ones(size(T_middle_phase1))* (min(peakMag.phase1) * 1.1);
                            plot(T_middle_phase1, yline_middle_phase1, 'LineWidth',7.5, 'Color','#83B4B3')

                            title("Manual segmentation Phase 1: " + fileName)

                           [loc, ~] = ginput(2);

                           if ~isempty(loc)
                               % find the selected locations that are close to
                               % the found peaks
                               ind(1) = interp1(peakLoc.phase1, 1:length(peakLoc.phase1) ,loc(1),'nearest');
                               ind(2) = interp1(peakLoc.phase1, 1:length(peakLoc.phase1), loc(2), 'nearest');

                               T_middle_phase1 = peakLoc.phase1(ind(1)):peakLoc.phase1(ind(2));
                               clear loc
                               close gcf


                               warning('Number of peaks found exceeds 7 \n 2 peaks manually selected phase 1')

                           else
                               error('Number of peaks found exceeds 7 \n NO peaks were selected phase 1')
                           end

                        end
                        %=========================================================

                        % throw error if the number of acceleration peaks does not reach 6
                        %=================================================================

                        if size(peakLoc.phase1, 1) < 6
                            fprintf('\t \t %s: Please select peaks phase 1 \n', content(file).name)

                            figure('Units','normalized','Position',[0 0 1 1]);
                            subplot(2,1,1);
                            % display the results of the change points

                            plot(df.pos.positionZ,"Color",[77 190 238]/255,"DisplayName","Input data")
                            hold on

                            % Plot segments between change points
                            plot(segmentMean,"Color",[64 64 64]/255,"DisplayName","Segment mean")

                            %Plot change points
                            x_rep = repelem(find(changeIndices),3);
                            y = repmat([ylim(gca) missing]',nnz(changeIndices),1);
                            plot(x_rep,y,"Color",[51 160 44]/255,"LineWidth",1,"DisplayName","Change points")
                            title("Number of change points: " + nnz(changeIndices))

                            xline(startPhase1, "Color", '#A2142F', "DisplayName",'StrPh1')
                            xline(endPhase1, "Color", '#A2142F', "DisplayName",'EndPh1')

                            xline(startPhase4, "Color", '#EDB120',"LineWidth",1, "DisplayName",'StrPh4')
                            xline(endPhase4, "Color", '#EDB120', "LineWidth",1,"DisplayName",'EndPh4')
                            hold off

                            %phase 1
                            subplot(2,1,2)
                            yline_phase1 = ones(size(T_phase1)) * (min(peakMag.phase1) * 1.2);

                            plot(df.SenAcc.SensorFreeX); hold on
                            plot(T_phase1, yline_phase1, 'LineWidth',10, 'Color', '#097392')
                            plot(grab_idx.phase1, peakMag.phase1(1:2:end),'o', 'MarkerSize', 7.5, 'MarkerFaceColor', 	'#A2142F', 'MarkerEdgeColor', 	'#A2142F')
                            plot(release_idx.phase1, peakMag.phase1(2:2:end), 'o', 'MarkerSize', 7.5, 'MarkerEdgeColor', 	'#0072BD', 'MarkerFaceColor', 	'#0072BD')

                            yline_middle_phase1 = ones(size(T_middle_phase1))* (min(peakMag.phase1) * 1.1);
                            plot(T_middle_phase1, yline_middle_phase1, 'LineWidth',7.5, 'Color','#83B4B3')

                            title("Manual segmentation Phase 1: " + fileName)

                           [loc, ~] = ginput(2);

                           if ~isempty(loc)
                               % find the selected locations that are close to
                               % the found peaks
                               ind(1) = interp1(peakLoc.phase1, 1:length(peakLoc.phase1) ,loc(1),'nearest');
                               ind(2) = interp1(peakLoc.phase1, 1:length(peakLoc.phase1), loc(2), 'nearest');

                               T_middle_phase1 = peakLoc.phase1(ind(1)):peakLoc.phase1(ind(2));
                               clear loc
                               close gcf


                               warning('Number of peaks found are less than 6 \n 2 peaks manually selected phase 1')

                           else
                               error('Number of peaks found are less than 6 \n NO peaks were selected phase 1')
                           end

                        end

                        %% middle phase 4
                        fprintf('\t \t %s: Segment middle arm movement phase 4 \n', content(file).name)

                        [peakLoc.phase4, peakMag.phase4] = peakfinder(df.SenAcc.SensorFreeX(T_phase4), [],[],-1);
                        peakLoc.phase4 = peakLoc.phase4 + startPhase4;
                        grab_idx.phase4 = peakLoc.phase4(1:2:end);
                        release_idx.phase4 = peakLoc.phase4(2:2:end);

                        T_middle_phase4 = grab_idx.phase4(2):grab_idx.phase4(3);

                        % throw error if the number of acceleration peaks
                        % exceeds 7
                        %================================================
                        if size(peakLoc.phase4, 1) >7
                            fprintf('\t \t %s: Please select peaks phase 4 \n', content(file).name)

                            figure('Units','normalized','Position',[0.1 0.1 0.75 0.75]);
                            subplot(2,1,1);
                            % display the results of the change points

                            plot(df.pos.positionZ,"Color",[77 190 238]/255,"DisplayName","Input data")
                            hold on

                            % Plot segments between change points
                            plot(segmentMean,"Color",[64 64 64]/255,"DisplayName","Segment mean")

                            %Plot change points
                            x_rep = repelem(find(changeIndices),3);
                            y = repmat([ylim(gca) missing]',nnz(changeIndices),1);
                            plot(x_rep,y,"Color",[51 160 44]/255,"LineWidth",1,"DisplayName","Change points")
                            title("Number of change points: " + nnz(changeIndices))

                            xline(startPhase1, "Color", '#A2142F', "DisplayName",'StrPh1')
                            xline(endPhase1, "Color", '#A2142F', "DisplayName",'EndPh1')

                            xline(startPhase4, "Color", '#EDB120',"LineWidth",1, "DisplayName",'StrPh4')
                            xline(endPhase4, "Color", '#EDB120', "LineWidth",1,"DisplayName",'EndPh4')
                            hold off

                            subplot(2,1,2)

                            plot(df.SenAcc.SensorFreeX); hold on

                            % phase 4
                            yline_phase4 = ones(size(T_phase4)) * (min(peakMag.phase4) * 1.2);
                            plot(T_phase4, yline_phase4, 'LineWidth',10, 'color', '#097392')

                            plot(grab_idx.phase4, peakMag.phase4(1:2:end),'o', 'MarkerSize', 7.5, 'MarkerFaceColor', 	'#A2142F', 'MarkerEdgeColor', 	'#A2142F')
                            plot(release_idx.phase4, peakMag.phase4(2:2:end), 'o', 'MarkerSize', 7.5, 'MarkerEdgeColor', 	'#0072BD', 'MarkerFaceColor', 	'#0072BD')

                            yline_middle_phase4 = ones(size(T_middle_phase4))*(min(peakMag.phase4) * 1.1);
                            plot(T_middle_phase4, yline_middle_phase4, 'LineWidth',7.5, 'Color','#83B4B3')

                            title("Manual segmentation Phase 4: " + fileName)

                            [loc, ~] = ginput(2);

                            if ~isempty(loc)
                                % find the selected locations that are close to
                                % the found peaks
                                ind(1) = interp1(peakLoc.phase4, 1:length(peakLoc.phase4) ,loc(1),'nearest');
                                ind(2) = interp1(peakLoc.phase4, 1:length(peakLoc.phase4), loc(2), 'nearest');

                                T_middle_phase4 = peakLoc.phase4(ind(1)):peakLoc.phase4(ind(2));
                                clear loc
                                close gcf


                                warning('Number of peaks found exceeds 7 \n peaks manually selected')
                            else
                                error('Number of peaks found exceeds 7 \n NO peaks were selected phase 1')
                            end   
                        end            
                        %================================================

                        % throw error if the number of acceleration peaks does not reach 6
                        %=================================================================

                        if size(peakLoc.phase4, 1) < 6 
                            fprintf('\t \t %s: Please select peaks phase 4 \n', content(file).name)

                            figure('Units','normalized','Position',[0.1 0.1 0.75 0.75]);
                            subplot(2,1,1);
                            % display the results of the change points

                            plot(df.pos.positionZ,"Color",[77 190 238]/255,"DisplayName","Input data")
                            hold on

                            % Plot segments between change points
                            plot(segmentMean,"Color",[64 64 64]/255,"DisplayName","Segment mean")

                            %Plot change points
                            x_rep = repelem(find(changeIndices),3);
                            y = repmat([ylim(gca) missing]',nnz(changeIndices),1);
                            plot(x_rep,y,"Color",[51 160 44]/255,"LineWidth",1,"DisplayName","Change points")
                            title("Number of change points: " + nnz(changeIndices))

                            xline(startPhase1, "Color", '#A2142F', "DisplayName",'StrPh1')
                            xline(endPhase1, "Color", '#A2142F', "DisplayName",'EndPh1')

                            xline(startPhase4, "Color", '#EDB120',"LineWidth",1, "DisplayName",'StrPh4')
                            xline(endPhase4, "Color", '#EDB120', "LineWidth",1,"DisplayName",'EndPh4')
                            hold off

                            subplot(2,1,2)
                            plot(df.SenAcc.SensorFreeX); hold on

                            % phase 4
                            yline_phase4 = ones(size(T_phase4)) * (min(peakMag.phase4) * 1.2);
                            plot(T_phase4, yline_phase4, 'LineWidth',10, 'color', '#097392')

                            plot(grab_idx.phase4, peakMag.phase4(1:2:end),'o', 'MarkerSize', 7.5, 'MarkerFaceColor', 	'#A2142F', 'MarkerEdgeColor', 	'#A2142F')
                            plot(release_idx.phase4, peakMag.phase4(2:2:end), 'o', 'MarkerSize', 7.5, 'MarkerEdgeColor', 	'#0072BD', 'MarkerFaceColor', 	'#0072BD')

                            yline_middle_phase4 = ones(size(T_middle_phase4))*(min(peakMag.phase4) * 1.1);
                            plot(T_middle_phase4, yline_middle_phase4, 'LineWidth',7.5, 'Color','#83B4B3')

                            title("Manual segmentation Phase 4: " + fileName)

                            [loc, ~] = ginput(2);

                            if ~isempty(loc)
                                % find the selected locations that are close to
                                % the found peaks
                                ind(1) = interp1(peakLoc.phase4, 1:length(peakLoc.phase4) ,loc(1),'nearest');
                                ind(2) = interp1(peakLoc.phase4, 1:length(peakLoc.phase4), loc(2), 'nearest');

                                T_middle_phase4 = peakLoc.phase4(ind(1)):peakLoc.phase4(ind(2));
                                clear loc
                                close gcf


                                warning('Number of peaks found are less than 6 \n peaks manually selected')
                            else
                                error('Number of peaks found are less than 6 \n NO peaks were selected phase 1')
                            end   
                        end

                        
                        %=================================================================
                        %% display the results
                        if plot_or_not
                            figure('Units','normalized','Position',[0 0 1 1]);
                            subplot(2,1,1);
                            % display the results of the change points

                            plot(df.pos.positionZ,"Color",[77 190 238]/255,"DisplayName","Input data")
                            hold on

                            % Plot segments between change points
                            plot(segmentMean,"Color",[64 64 64]/255,"DisplayName","Segment mean")

                            %Plot change points
                            x_rep = repelem(find(changeIndices),3);
                            y = repmat([ylim(gca) missing]',nnz(changeIndices),1);
                            plot(x_rep,y,"Color",[51 160 44]/255,"LineWidth",1,"DisplayName","Change points")
                            title("Number of change points: " + nnz(changeIndices))

                            xline(startPhase1, "Color", '#A2142F', "DisplayName",'StrPh1')
                            xline(endPhase1, "Color", '#A2142F', "DisplayName",'EndPh1')

                            xline(startPhase4, "Color", '#EDB120',"LineWidth",1, "DisplayName",'StrPh4')
                            xline(endPhase4, "Color", '#EDB120', "LineWidth",1,"DisplayName",'EndPh4')
                            hold off

                            %phase 1
                            subplot(2,1,2)
                            yline_phase1 = ones(size(T_phase1)) * (min(peakMag.phase1) * 1.2);

                            plot(df.SenAcc.SensorFreeX); hold on
                            plot(T_phase1, yline_phase1, 'LineWidth',10, 'Color', '#097392')
                            plot(grab_idx.phase1, peakMag.phase1(1:2:end),'o', 'MarkerSize', 7.5, 'MarkerFaceColor', 	'#A2142F', 'MarkerEdgeColor', 	'#A2142F')
                            plot(release_idx.phase1, peakMag.phase1(2:2:end), 'o', 'MarkerSize', 7.5, 'MarkerEdgeColor', 	'#0072BD', 'MarkerFaceColor', 	'#0072BD')

                            yline_middle_phase1 = ones(size(T_middle_phase1))* (min(peakMag.phase1) * 1.1);
                            plot(T_middle_phase1, yline_middle_phase1, 'LineWidth',7.5, 'Color','#83B4B3')

                            % phase 4
                            yline_phase4 = ones(size(T_phase4)) * (min(peakMag.phase4) * 1.2);
                            plot(T_phase4, yline_phase4, 'LineWidth',10, 'color', '#097392')

                            plot(grab_idx.phase4, peakMag.phase4(1:2:end),'o', 'MarkerSize', 7.5, 'MarkerFaceColor', 	'#A2142F', 'MarkerEdgeColor', 	'#A2142F')
                            plot(release_idx.phase4, peakMag.phase4(2:2:end), 'o', 'MarkerSize', 7.5, 'MarkerEdgeColor', 	'#0072BD', 'MarkerFaceColor', 	'#0072BD')

                            yline_middle_phase4 = ones(size(T_middle_phase4))*(min(peakMag.phase4) * 1.1);
                            plot(T_middle_phase4, yline_middle_phase4, 'LineWidth',7.5, 'Color','#83B4B3')

                            title("visual segmentation: " + fileName)

                        end
                        disp('      ')

                        %% extract relevant kinematics
                        fprintf('\t\t %s: Extract relevant kinematics \n', fileName)
                        % initialise joint names
                        jointNames = ['Scapula', "Glenohumeraal", "Elbow", "Wrist"];
                        for jnt = 1:4
                            % initiation of joint names
                            IK_X = [jointNames{jnt}, '_abbuction'];
                            IK_Y = [jointNames{jnt}, '_rotation'];
                            IK_Z = [jointNames{jnt}, '_flexion'];

                            Data_out.(Timepoint).(subj_name).phase1.(arm).(IK_X)(:, counter) = normalisation(jointData(jointNo(jnt)).jointAngle(:, 1), T_middle_phase1);
                            Data_out.(Timepoint).(subj_name).phase1.(arm).(IK_Y)(:, counter) = normalisation(jointData(jointNo(jnt)).jointAngle(:, 2), T_middle_phase1);
                            Data_out.(Timepoint).(subj_name).phase1.(arm).(IK_Z)(:, counter) = normalisation(jointData(jointNo(jnt)).jointAngle(:, 3), T_middle_phase1);


                            Data_out.(Timepoint).(subj_name).phase4.(arm).(IK_X)(:, counter) = normalisation(jointData(jointNo(jnt)).jointAngle(:, 1), T_middle_phase4);
                            Data_out.(Timepoint).(subj_name).phase4.(arm).(IK_Y)(:, counter) = normalisation(jointData(jointNo(jnt)).jointAngle(:, 2), T_middle_phase4);
                            Data_out.(Timepoint).(subj_name).phase4.(arm).(IK_Z)(:, counter) = normalisation(jointData(jointNo(jnt)).jointAngle(:, 3), T_middle_phase4);

                        end
                        % Trunk data
                        trunk = jointData(2).jointAngle + jointData(3).jointAngle + jointData(4).jointAngle; %+ jointData(5).jointAngle;
                        Trunk_lateroflexion = trunk(:,1);
                        Trunk_rotation = trunk(:,2);
                        Trunk_flexion = trunk(:,3);


                        % initiation of trunk names
                        IK_X = 'Trunk_lateroflexion';
                        IK_Y = 'Trunk_rotation';
                        IK_Z = 'Trunk_flexion';

                        Data_out.(Timepoint).(subj_name).phase1.(arm).(IK_X)(:, counter) = normalisation(trunk(:,1), T_middle_phase1);
                        Data_out.(Timepoint).(subj_name).phase1.(arm).(IK_Y)(:, counter) = normalisation(trunk(:, 2), T_middle_phase1);
                        Data_out.(Timepoint).(subj_name).phase1.(arm).(IK_Z)(:, counter) = normalisation(trunk(:, 3), T_middle_phase1);


                        Data_out.(Timepoint).(subj_name).phase4.(arm).(IK_X)(:, counter) = normalisation(jointData(jointNo(jnt)).jointAngle(:, 1), T_middle_phase4);
                        Data_out.(Timepoint).(subj_name).phase4.(arm).(IK_Y)(:, counter) = normalisation(jointData(jointNo(jnt)).jointAngle(:, 2), T_middle_phase4);
                        Data_out.(Timepoint).(subj_name).phase4.(arm).(IK_Z)(:, counter) = normalisation(jointData(jointNo(jnt)).jointAngle(:, 3), T_middle_phase4);
                       
                        %% clear certain variables
                        clear peakLoc peakMag

                    end% first check
                end % end if movement && .mvnx
            catch ME
                % this section of code displays the subject and filenames
                % on which an error occured
                %                 formatSpec = '\t ERROR: %s at TIMEPOINT: %s in FILENAME: %s \n %s LINE: %s ';
                %                 fprintf(formatSpec, subj_name, Timepoint, fileName, ME.message, num2str(ME.stack.line));
                %                 disp(' ')

                cprintf('*text', '\t\t ERROR: %s at TIMEPOINT: %s in FILENAME: %s \n \t\t %s LINE: %s \n', ...
                    subj_name, Timepoint, fileName, ME.message, num2str(ME.stack.line))

            end
        end %end number of files

        %this section of code writes the peaks information to an excel file per subject
        % This also indicates whether the ULIFT task is segmented (1) or not
        % (0)

        %         if strcmp(Timepoint, 'T0')
        %             T = table(ppID, filename, Tpoint, Phase1, Phase4, run);
        %             writetable(T, fullfile(path.table, [subj_name, '.xlsx']))
        %         else
        %             tablefile = fullfile(path.table, [subj_name, '.xlsx']);
        %             T_orig = readtable(tablefile);
        %
        %             T_new = table(ppID, filename, Tpoint, Phase1, Phase4, run);
        %             T = [T_orig; T_new];
        %             writetable(T, fullfile(path.table, [subj_name, '.xlsx']))
        %         end
    end% end check if subject path exists

    %% Save data
    %--------------

    %     if exist('Data_out','var')
    %         if exist(path.out,'file')
    %             load(path.out)
    %         end
    %
    %         [Data.(subj_name).ULIFT.(Timepoint)] = Data_out.ULIFT.(Timepoint);
    %         save(path.out,'Data')
    %         clear Data Data_out
    %     end

    fprintf('*********Finished %s **********\n', subj_name)

end %end subjects
cd("C:\Users\u0117545\Documents\GitHub\ULIFT_BC")

%% plot joint angle data
JointAngles = fieldnames(Data_out.(Timepoint).(subj_name).phase4.(arm));

figure
p = tiledlayout('flow', 'TileSpacing','compact');
title(p,['Phase 1 Average left and right ', (Timepoint)])
xlabel(p,'% movement')
ylabel(p,['Joint angle', char(176)])

for jnt = 1:length(JointAngles)
    L.Phase1.(JointAngles{jnt}) = Data_out.(Timepoint).(subj_name).phase4.L.(JointAngles{jnt});
    R.Phase1.(JointAngles{jnt}) = Data_out.(Timepoint).(subj_name).phase4.R.(JointAngles{jnt});

    nexttile
    plot(L.Phase1.(JointAngles{jnt}), 'DisplayName', "Avg Left")
    hold on
    plot(R.Phase1.(JointAngles{jnt}), 'DisplayName', "Avg Right")
    xlim([0, 100])

    [PlotTitle] = regexprep(JointAngles{jnt}, '_', ' ');
    title(PlotTitle)
end


diary off
