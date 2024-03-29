




clearvars; close all; clc
%% Input data paths and names
% 
%  INPUT 1: corresponding path where all files are to path.root
%  INPUT 2: corresponging ppID number that is being analysed 
%  INPUT 3: approximate timerange when this participant was measured
%  INPUT 4: file names to the data for left, right, and video files respectively

path.root   = 'C:\Users\u0117545\KU Leuven\An De Groef - DATA\LUM'; 
ppID        = 'L_001'; 

subj_path   = fullfile(path.root, ppID, 'csv');

content     = dir(subj_path);
nFiles      = length(content);

for file = 1:nFiles
    %disp(content(file).name)

    if contains(content(file).name, 'L.2')
        disp(content(file).name)
        disp("  File for the LEFT arm")

        file_l = content(file).name;
    elseif contains(content(file).name, 'R.2')
        disp(content(file).name)
        disp("  File for the RIGHT arm")

        file_r = content(file).name;
    end

end

file_l = 'MOS2D31170588 (2022-11-18).csv';
file_r = 'MOS2E08190543 (2022-11-18).csv';

%% manulal input no longer nesesary after this point
fileName_L      = fullfile(subj_path, file_l);
fileName_R      = fullfile(subj_path, file_r);

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 4);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Var1", "X", "Y", "Z"];
opts.SelectedVariableNames = ["X", "Y", "Z"];
opts.VariableTypes = ["string", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "Var1", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Var1", "EmptyFieldRule", "auto");

R = readtable(fileName_R, opts);
L = readtable(fileName_L, opts);

% Clear temporary variables
clear opts


%% calculate timedelay between the signals

% fig = figure("Units","normalized", 'outerposition', [0,0,1,1]);
% sph = stackedplot(R);
% 
% % Assign callback function when user clicks on figure
% fig.WindowButtonDownFcn = {@myTR, sph};

% 


if strcmp(ppID, 'L_009')
    figure("Units","normalized", 'outerposition', [0,0,1,1]);
    subplot(3,1,1); plot(L.X(height(L)/2:end))
    subplot(3,1,2); plot(L.Y(height(L)/2:end))
    subplot(3,1,3); plot(L.Z(height(L)/2:end))
    [gi,~] = ginput(2);

    gi = gi+ (height(L)/2);
else
    figure("Units","normalized", 'outerposition', [0,0,1,1]);
    subplot(3,1,1); plot(L.X(1:height(L)/2))
    subplot(3,1,2); plot(L.Y(1:height(L)/2))
    subplot(3,1,3); plot(L.Z(1:height(L)/2))
    [gi,~] = ginput(2);
end


signal1 = L.X(round(gi(1)):round(gi(2)));
signal2 = R.X(round(gi(1)):round(gi(2)));

corr = xcorr(signal1, signal2);
[~, idx] = max(corr);
delay = idx - length(signal1);

signal2_synced = circshift(signal2, delay);

R_synced = circshift(R, delay);


figure; 
plot(L.X(round(gi(1)):round(gi(2))));
hold on
plot(R.X(round(gi(1)):round(gi(2))));
plot(R_synced.X(round(gi(1)):round(gi(2))))
%% Create timetable data and cut out relevant time section
% Determine timerange based on right arm.

% L_time = table2timetable(L);

fig = figure("Units","normalized", 'outerposition', [0,0,1,1]);
sph = stackedplot(L);
title("Left file")

% Assign callback function when user clicks on figure
fig.WindowButtonDownFcn = {@myTR, sph};


%R_time = table2timetable(R);

fig = figure("Units","normalized", 'outerposition', [0,0,1,1]);
sph = stackedplot(R_synced);
title('Right file')
% Assign callback function when user clicks on figure
fig.WindowButtonDownFcn = {@myTR, sph};





%% Footnotes
% [1] The pause(0.5) is needed to give time for the data tip values to 
%   update. Notice that the line must remain still for a brief period
%   before the datatips appear.  In that time, the datatip values from
%   the previous datatip appearance are stored and we must wait for them
%   to regenerate before returning the updated value.  This lag is a 
%   property of the stackedplot data cursor:  Sdc.LingerTime. 

function myTR(~,~,sph)
% sph is the stackedplot handle.
pause(0.5) % See footnote [1]
% Avoid seeing warning associated with undocumented methods
origState = warning('query', 'MATLAB:structOnObject');
cleanup = onCleanup(@()warning(origState));
warning('off','MATLAB:structOnObject')
S = struct(sph);            % Undocumented
Sdc = struct(S.DataCursor); % Undocumented
% Get data tip values
dataTipLabels = Sdc.DatatipLabels;
dataTipStr = get([dataTipLabels{:}],'String');
axLabels = sph.DisplayLabels;

% If the vertical CursorLine is not visible (ie, the mouse is not in the axes) or
% if the data tips are all NaN (sometimes this happens on the first click only)
% then do nothing.  Otherwise, print the data tip values and axis labels to command
% window.
if strcmpi(Sdc.CursorLine.Visible,'off') || all(cellfun('isempty',dataTipStr(1:end-1)))
    return
else
    yDataStr = compose('   %s = %.3g', string(axLabels), string(dataTipStr(1:end-1)));
    fprintf('\nAt x = %s\n%s\n',dataTipStr{end},strjoin(yDataStr,newline))
end
end


