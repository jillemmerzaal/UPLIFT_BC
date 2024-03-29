%% timepoints of interst

timepoints = readtable("C:\Users\u0117545\Documents\GitHub\ULIFT_BC\ValidationStartEnd.xlsx");


timepoints = timepoints(strcmp(timepoints.subj_id, subj_name),:);

interest = timepoints(strcmp(timepoints.Trial, fileName),:);

x1 = interest.Start_phase1;
x2 = interest.Stop_phase1;
x3 = interest.Start_phase4;
x4 = interest.Stop_phase4;

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
% freqz(b,a)

position = filtfilt(b,a, segmentData(segmentno).position);
positionX = position(:,1);
positionY = position(:,2);
positionZ = position(:,3);
positionVec = vecnorm(position, 2,2);


velocity = filtfilt(b,a, segmentData(segmentno).velocity);
velocityX = velocity(:,1);
velocityY = velocity(:,2);
velocityZ = velocity(:,3);
velocityVec = vecnorm(velocity, 2, 2);

acceleration = filtfilt(b,a, segmentData(segmentno).acceleration);
accelerationX = acceleration(:,1);
accelerationY = acceleration(:,2);
accelerationZ = acceleration(:,3);
accelerationVec = vecnorm(acceleration, 2, 2);

angularVel_LA = filtfilt(b,a, segmentData(segmentno).angularVelocity);
angularVelX_LA = angularVel_LA(:,1);
angularVelY_LA = angularVel_LA(:,2);
angularVelZ_LA = angularVel_LA(:,3);
angularVelVec = vecnorm(angularVel_LA, 2, 2);
angularVelDiff = [diff(angularVelVec); 0];

% angularVelX_UA = segmentData(12).angularVelocity(:,1);
% angularVelY_UA = segmentData(12).angularVelocity(:,2);
% angularVelZ_UA = segmentData(12).angularVelocity(:,3);

SensorFree = filtfilt(b,a, sensorData(sensorno).sensorFreeAcceleration);
SensorFreeX = SensorFree(:,1);
SensorFreeY = SensorFree(:,2);
SensorFreeZ = SensorFree(:,3);
SensorFreeVec = vecnorm(SensorFree,2,2);
SensorFreeDiff = [diff(SensorFreeVec); 0];
%SensorFreeDiff = [SensorFreeDiff; 0];
% figure; plot(tp)

%% dataframe
%df.acc      = table(accelerationX, accelerationY, accelerationZ, accelerationVec);
%df.vel      = table(velocityX, velocityY, velocityZ, velocityVec);
%df.pos      = table(positionX, positionY, positionZ, positionVec);
df.Avel     = table(angularVelX_LA, angularVelY_LA, angularVelZ_LA, angularVelVec, angularVelDiff);
df.SenAcc   = table(SensorFreeX, SensorFreeY, SensorFreeZ, SensorFreeVec, SensorFreeDiff);

%% plot full signals.


%
% for nPlot = 1:length(signals)
%     figure;
%     h = stackedplot(df.(signals{nPlot}));
%     title(plottitle(nPlot))
%
%     %based on xsens
%     ax = findobj(h.NodeChildren, 'Type','Axes');
%     arrayfun(@(h)xline(h,x1,'LineWidth',1.5, "Color", '#A2142F', "DisplayName",'StrPh1'),ax)
%     arrayfun(@(h)xline(h,x2, 'LineWidth', 1.5,"Color", '#A2142F', "DisplayName",'EndPh1'), ax)
%     arrayfun(@(h)xline(h,x3,'LineWidth',1.5, "Color", '#A2142F', "DisplayName",'StrPh4'),ax)
%     arrayfun(@(h)xline(h,x4, 'LineWidth', 1.5,"Color", '#A2142F', "DisplayName",'EndPh4'), ax)
%
%     %based on zerocrossing velocity Z
%     arrayfun(@(h)xline(h,startPhase1,'LineWidth',1.5,'LineStyle', ':', "Color", 'green', "DisplayName",'StrPh1'),ax)
%     arrayfun(@(h)xline(h,endPhase1, 'LineWidth', 1.5, 'LineStyle', ':', "Color", 'green', "DisplayName",'EndPh1'), ax)
%     arrayfun(@(h)xline(h,startPhase4,'LineWidth',1.5,'LineStyle', ':', "Color", 'green', "DisplayName",'StrPh4'),ax)
%     arrayfun(@(h)xline(h,endPhase4, 'LineWidth', 1.5, 'LineStyle', ':', "Color", 'green', "DisplayName",'EndPh4'), ax)
%
%     % change indices
%     arrayfun(@(h)xline(h,x(1),'LineWidth',1.5,'LineStyle', ':', "Color", 'black', 'DisplayName','change1'),ax)
%     arrayfun(@(h)xline(h,x(2), 'LineWidth', 1.5, 'LineStyle', ':', "Color", 'black', 'DisplayName','change2'), ax)
%
% end

%% find end first phase
% start = 1;
% stop = x(1);
%
% for nPlot = 1:length(signals)
%     figure;
%     h = stackedplot(df.(signals{nPlot})(start:stop,:));
%     title(plottitle(nPlot))
%
%
%     %based on xsens
%     ax = findobj(h.NodeChildren, 'Type','Axes');
%     arrayfun(@(h)xline(h,x1,'LineWidth',1.5, "Color", '#A2142F', "DisplayName",'StrPh1'),ax)
%     arrayfun(@(h)xline(h,x2, 'LineWidth', 1.5,"Color", '#A2142F', "DisplayName",'EndPh1'), ax)
%
%
%     %based on zerocrossing velocity Z
%     arrayfun(@(h)xline(h,startPhase1,'LineWidth',1.5,'LineStyle', ':', "Color", 'green', "DisplayName",'StrPh1'),ax)
%     arrayfun(@(h)xline(h,endPhase1, 'LineWidth', 1.5, 'LineStyle', ':', "Color", 'green', "DisplayName",'EndPh1'), ax)
%
%
%     % change indices
%     arrayfun(@(h)xline(h,x(1),'LineWidth',1.5,'LineStyle', ':', "Color", 'black', 'DisplayName','change1'),ax)
% end

%% New start phase 1
[temp, P] = islocalmax(df.SenAcc.SensorFreeDiff(1:x(1)));


localmax.all = temp;
stderror= std(P(localmax.all));
average = mean(P(localmax.all));
%Thresh = mean(P(localmax.all)) - mean(P(localmax.all))*0.25;
Thresh = average*0.25;

clear temp P
[temp, P] = islocalmax(df.SenAcc.SensorFreeDiff(1:x(1)), 'MinProminence',Thresh);
localmax.thresh = temp;

N = 1:height(df.SenAcc);

figure;
plot(N,df.SenAcc.SensorFreeDiff,N(localmax.all),df.SenAcc.SensorFreeDiff(localmax.all),'r*')
hold on
xline(x(1))
plot(N,df.SenAcc.SensorFreeDiff,N(localmax.thresh),df.SenAcc.SensorFreeDiff(localmax.thresh),'g*')

localmax.prominent = N(localmax.thresh)
localmax.incon = N(localmax.all)

startPhase1_new = localmax.prominent(1);

clear localmax temp P

%% New start phase 1 version 2 
% [temp, P] = islocalmax(df.Avel.angularVelDiff(1:x(1)));
% localmax.all = temp;
% average = mean(P(localmax.all));
% Thresh = average *1.5;


[peakLoc, peakMag] = peakfinder(df.Avel.angularVelDiff(1:x(1)));
localmax.all = peakLoc;
Thresh = mean(peakMag) *1.5;
[peakLoc, peakMag] = peakfinder(df.Avel.angularVelDiff(1:x(1)), [], Thresh);
localmax.thresh = peakLoc;



% clear temp P
% [temp, P] = islocalmax(df.Avel.angularVelDiff(1:x(1)), 'MinProminence',Thresh);
% localmax.thresh = temp;

N = 1:height(df.Avel);

figure;
plot(N,df.Avel.angularVelDiff,localmax.all,df.Avel.angularVelDiff(localmax.all),'r*')
hold on
xline(x(1))
% plot(N,df.Avel.angularVelDiff,N(localmax.thresh),df.Avel.angularVelDiff(localmax.thresh),'g*')

% localmax.prominent = N(localmax.thresh)
% localmax.incon = N(localmax.all)

if isempty(localmax.thresh)
    startPhase1_new2 = localmax.all(1);
else
    startPhase1_new2 = localmax.thresh(1);
end
clear localmax temp P

%% New end phase 1
[temp, P] = islocalmin(df.SenAcc.SensorFreeX(1:x(1)));

localmin.all = temp;
Thresh = mean(P(localmin.all));

clear temp P
[temp, P] = islocalmin(df.SenAcc.SensorFreeX(1:x(1)), 'MaxNumExtrema', 12, 'MinProminence',Thresh);
localmin.thresh = temp
N = 1:height(df.SenAcc);

figure;
plot(N,df.SenAcc.SensorFreeX,N(localmin.all),df.SenAcc.SensorFreeX(localmin.all),'r*')
hold on
xline(x(1))

%select the less prominent minima between the last two most prominent minima
localmin.prominent = N(localmin.thresh)
localmin.incon = N(localmin.all)

% endPhase1_new = localmin.incon(find(localmin.incon == localmin.prominent(end))-1)
endPhase1_new = localmin.incon(end-1);

% idx = N(localmin);
%
% grab_phase1 = idx(1:2:end);
% go_phase1 = idx(2:2:end);
% endPhase1_new = go_phase1(end);
%
% clear localmin P localmax P_max idx
clear localmin P temp


%% new start phase 4
temp = df.SenAcc.SensorFreeX(x(2):end);
[min, P] = islocalmin(temp);

localmin.all = min;
Thresh = mean(P(localmin.all)) + std(P(localmin.all)) *0.1;
clear min P
[min, P] = islocalmin(temp, 'MinProminence',Thresh);
localmin.thresh = min;


figure;
plot(N,df.SenAcc.SensorFreeX,N(localmin.all)+x(2),temp(localmin.all),'r*')

%
localmin.prominent = N(localmin.thresh)+ x(2);
localmin.incon = N(localmin.all)+ x(2);
startPhase4_new = localmin.prominent(1);

% idx = N(localmin) + x(2);
% go_phase4 = idx(1:2:end);
% grab_phase4 = idx(2:2:end);

% startPhase4_new = grab_phase4(1);

clear localmin min P

%% new end phase 4
% temp = df.SenAcc.SensorFreeDiff(x(2):end);
% [min, P] = islocalmin(temp);
%
% localmin.all = min;
% Thresh = mean(P(localmin.all)) + std(P(localmin.all)) *0.1;
% clear min P
% [min, P] = islocalmin(temp, 'MinProminence',Thresh);
% localmin.thresh = min
%
% figure;
% plot(N,df.SenAcc.SensorFreeDiff,N(localmin.all)+x(2),temp(localmin.all),'r*')
%
% localmin.prominent = N(localmin.thresh)+ x(2)
% localmin.incon = N(localmin.all)+ x(2)
%
% endPhase4_new = localmin.prominent(end)
% clear localmin temp P min
%% new end phase 4 version 2
temp =df.SenAcc.SensorFreeY(x(2):end);
[min, P] = islocalmin(temp);

localmin.all = min;
Thresh = mean(P(localmin.all)) + std(P(localmin.all)) *0.1;
clear min P

[min, P] = islocalmin(temp, 'MinProminence',Thresh);
localmin.thresh = min;
clear min P

figure;
plot(N, df.SenAcc.SensorFreeY, N(localmin.all)+x(2), temp(localmin.all), 'r*')

localmin.prominent = N(localmin.thresh)+x(2);
localmin.incon = N(localmin.all) + x(2);

endPhase4_new2 = localmin.prominent(end);

%% new end phase 4 version 3
temp =df.Avel.angularVelDiff(x(2):end)*-1;
[peakLoc, peakMag] = peakfinder(temp);
localmin.all = peakLoc + x(2);

% localmin.all = min;
average = mean(peakMag);
Thresh = average *1.5;

% Thresh = mean(P(localmin.all)) + std(P(localmin.all)) *0.25;
% clear min P

% [min, P] = islocalmin(temp, 'MinProminence',Thresh);
[peakLoc, peakMag] = peakfinder(temp, [], Thresh);


localmin.thresh = peakLoc + x(2);

% localmin.thresh = min;
% clear min P

% figure;
% plot(N, df.Avel.angularVelDiff, N(localmin.all)+x(2), temp(localmin.all), 'r*')
% hold on
% plot(N, df.Avel.angularVelDiff, N(localmin.thresh)+x(2), temp(localmin.thresh), 'g*')


figure;
plot(N,df.Avel.angularVelDiff,localmin.thresh,df.Avel.angularVelDiff(localmin.thresh),'r*')
hold on
xline(x(1))



localmin.prominent = N(localmin.thresh)+x(2);
localmin.incon = N(localmin.all) + x(2);

% endPhase4_new3 = localmin.prominent(end);
endPhase4_new3 = peakLoc(end) + x(2);
%% check the shit

signals = fieldnames(df);
for nPlot = 1:length(signals)
    plottitle = {[signals{nPlot} ' data ' fileName]};

    figure;
    h = stackedplot(df.(signals{nPlot}));
    title(plottitle)

    %based on xsens
    ax = findobj(h.NodeChildren, 'Type','Axes');
    arrayfun(@(h)xline(h,x1,'LineWidth',1.5, "Color", '#A2142F', "DisplayName",'StrPh1'),ax)
    arrayfun(@(h)xline(h,x2, 'LineWidth', 1.5,"Color", '#A2142F', "DisplayName",'EndPh1'), ax)
    arrayfun(@(h)xline(h,x3,'LineWidth',1.5, "Color", '#A2142F', "DisplayName",'StrPh4'),ax)
    arrayfun(@(h)xline(h,x4, 'LineWidth', 1.5,"Color", '#A2142F', "DisplayName",'EndPh4'), ax)

    %based on zerocrossing velocity Z
    arrayfun(@(h)xline(h,startPhase1,'LineWidth',1.5, "Color", 'green', "DisplayName",'StrPh1'),ax)
    arrayfun(@(h)xline(h,endPhase1, 'LineWidth', 1.5,  "Color", 'green', "DisplayName",'EndPh1'), ax)
    arrayfun(@(h)xline(h,startPhase4,'LineWidth',1.5, "Color", 'green', "DisplayName",'StrPh4'),ax)
    arrayfun(@(h)xline(h,endPhase4, 'LineWidth', 1.5,  "Color", 'green', "DisplayName",'EndPh4'), ax)

    %based on sensor free acceleration in X direction
    arrayfun(@(h)xline(h,endPhase1_new, 'LineWidth', 1.5, 'LineStyle', ':', "Color", 'Blue', "DisplayName",'EndPh1'), ax)
    arrayfun(@(h)xline(h, startPhase4_new, 'LineWidth', 1.5, 'LineStyle', ':', 'Color', 'Blue', 'DisplayName', 'StrPh4'), ax)

    %new start phase 1 based on the highest rate of change of the
    %acceleration vector
    arrayfun(@(h)xline(h, startPhase1_new, 'Linewidth', 1.5, 'LineStyle', ':', 'Color', 'Blue', 'DisplayName', 'StartPh1'), ax)
    arrayfun(@(h)xline(h, endPhase4_new2, 'Linewidth', 1.5, 'LineStyle', ':', 'Color', 'Blue', 'DisplayName', 'endPh4'), ax)

    %new version 2 start phase 1 based on the highest rate of change of the
    %angular velocity vector
    arrayfun(@(h)xline(h, startPhase1_new2, 'Linewidth', 1.5, 'LineStyle', ':', 'Color', 'm', 'DisplayName', 'StartPh1'), ax)

    arrayfun(@(h)xline(h, endPhase4_new3, 'Linewidth', 1.5, 'LineStyle', ':', 'Color', 'm', 'DisplayName', 'StartPh1'), ax)

    % change indices
    arrayfun(@(h)xline(h,x(1),'LineWidth',1.5,'LineStyle', ':', "Color", 'black', 'DisplayName','change1'),ax)
    arrayfun(@(h)xline(h,x(2), 'LineWidth', 1.5, 'LineStyle', ':', "Color", 'black', 'DisplayName','change2'), ax)

  
end


%% plot for paper? 

plottable = table(df.pos.positionZ, df.Avel.angularVelDiff, df.SenAcc.SensorFreeX);
plottable.Properties.VariableNames = {'Position Z', 'Diff Angular Vel Vec', 'Acc X'};
writetable(plottable, 'EventDetectionTimecures.csv')

interest.Start_phase1(2)    = startPhase1_new2;
interest.Stop_phase1(2)     = endPhase1_new;
interest.Start_phase4(2)    = startPhase4_new;
interest.Stop_phase4(2)     = endPhase4_new3;
interest.Change1(2)         = find(changeIndices, 1, 'first');
interest.Change2(2)         = find(changeIndices, 1, 'last');


writetable(interest, 'EventDetectionPOI.csv')

