%% SPM preperation ULIFT data
%
% This is the second step in the analysis of the ULIFT data. This code
% works completely automatically. 
% Fill out the path where the code and the data is stored.
% Fill out the timepoints and phase you want to analyse. The timepoints can
% be multiple at once. But just fill out 1 phase.
% Set the number of participants. 


%% FIX THIS: Jill 
% 
% at the end of the script, make the table also in a loop. 


clearvars; clc; close all

%% set input path
% manual input needed for paths
path.code   = "C:\Users\u0117545\Documents\GitHub\ULIFT_BC\ULIFT";
path.root   = 'C:\Users\u0117545\KU Leuven\An De Groef - DATA';

% no manual input needed
cd(path.code)
addpath(path.code)
path.in     = fullfile(path.code, 'Output', 'Database_ULIFT.mat');
path.out    = fullfile(path.code,'Output','Database_SPM.mat');


%% load data
if exist(path.in,'file')
    load(path.in)
end
affected_table = readtable(fullfile(path.root, "Aangedane zijde.xlsx"));

%% set up
Timepoints = {'T0', 'T1'};
Phase       = 'phase4';


jointsOfInterst ={'Scapula_abbuction'
    'Scapula_rotation'
    'Scapula_flexion'
    'Glenohumeraal_abbuction'
    'Glenohumeraal_rotation'
    'Glenohumeraal_flexion'
    'Elbow_flexion'
    'Trunk_lateroflexion'
    'Trunk_rotation'
    'Trunk_flexion'};

nangles = size(jointsOfInterst,1);
ntime = size(Timepoints,2);


for t = 1:ntime
    for subj = (1:20) 
        if subj < 10
            subj_name   = ['BC_00' num2str(subj)];
        elseif subj < 100
            subj_name   = ['BC_0' num2str(subj)];
        else
            subj_name   = ['BC_', num2str(subj)];
        end

        %% run data
        % check if participant and timepoint are in struct
        if isfield(Data, subj_name)
            if isfield(Data.(subj_name), Timepoints{t})
                % find the affected side
                idx = find(strcmp(affected_table.ppID, subj_name));
                involved = affected_table(idx,:).involved;

                % setup row number for exporting the data
                d = strfind(subj_name,'_');
                rownr = str2double(subj_name(d+1:end));
%                 fprintf('\n')
%                 cprintf('*text', 'Processing: %s at Timepoint: %s....... \n', subj_name, Timepoints{t})
%                 cprintf('*text', '\t\t rownumber: %d \n', rownr)

                %% Average joint data from the joints of interest
                for ang = 1:nangles
                    % if RIGHT is affected
                    if strcmp(involved, 'R')
                        if strcmp(Phase, 'phase1')
                            if strcmp(jointsOfInterst{ang}, 'Trunk_lateroflexion') || strcmp(jointsOfInterst{ang}, 'Trunk_rotation')
                                temp.df = mean(-1 .* Data.(subj_name).(Timepoints{t}).(Phase).R.(jointsOfInterst{ang}), 2)';
                                BC.aff.(Phase).(Timepoints{t}).(jointsOfInterst{ang})(rownr,:) = table(string(subj_name), string(Timepoints{t}), temp.df);
                                clear temp
                            else
                                temp.df = mean(Data.(subj_name).(Timepoints{t}).(Phase).R.(jointsOfInterst{ang}), 2)';
                                BC.aff.(Phase).(Timepoints{t}).(jointsOfInterst{ang})(rownr,:) = table(string(subj_name), string(Timepoints{t}), temp.df);
                                clear temp
                            end
                        elseif strcmp(Phase, 'phase4')
                            temp.df = mean(Data.(subj_name).(Timepoints{t}).(Phase).R.(jointsOfInterst{ang}), 2)';
                            BC.aff.(Phase).(Timepoints{t}).(jointsOfInterst{ang})(rownr,:) = table(string(subj_name), string(Timepoints{t}), temp.df);
                            clear temp
                        end
                        % BC.unaff.(Phase).(Timepoints{t}).(jointsOfInterst{ang})(rownr,:) = mean(Data.(subj_name).(Timepoints{t}).(Phase).L.(jointsOfInterst{ang}), 2)';


                        % if LEFT is affected
                    elseif strcmp(involved, 'L')
                        temp.df =  mean(Data.(subj_name).(Timepoints{t}).(Phase).L.(jointsOfInterst{ang}), 2)';
                        BC.aff.(Phase).((Timepoints{t})).(jointsOfInterst{ang})(rownr,:) = table(string(subj_name), string(Timepoints{t}), temp.df);

                        %                        if strcmp(Phase, 'phase1')
                        %                            if strcmp(jointsOfInterst{ang}, 'Trunk_lateroflexion') || strcmp(jointsOfInterst{ang}, 'Trunk_rotation')
                        %                                BC.unaff.(Phase).(Timepoints{t}).(jointsOfInterst{ang})(rownr,:) = mean(-1 .* Data.(subj_name).(Timepoints{t}).(Phase).R.(jointsOfInterst{ang}), 2)';
                        %                            else
                        %                                BC.unaff.(Phase).(Timepoints{t}).(jointsOfInterst{ang})(rownr,:) = mean(Data.(subj_name).(Timepoints{t}).(Phase).R.(jointsOfInterst{ang}), 2)';
                        %                            end
                        %                        elseif strcmp(Phase, 'phase4')
                        %                            BC.unaff.(Phase).(Timepoints{t}).(jointsOfInterst{ang})(rownr,:) = mean(Data.(subj_name).(Timepoints{t}).(Phase).R.(jointsOfInterst{ang}), 2)';
                        %                        end
                        clear temp

                    end
                end

               
            end % check if timepoint is in that particular subject.
        end % check if subject name is in struct
    end % number of subjects

    %% write data to excel
    for ang = 1:nangles
        temp.table = [BC.aff.(Phase).(Timepoints{t}).(jointsOfInterst{ang})];
        temp.table.Properties.VariableNames = {'ppID', 'Time', jointsOfInterst{ang}};

        % create filename
        fileName = fullfile(path.code, 'Output', [Timepoints{t}, '_', Phase, '.xlsx']);

        %  write the table to an excel file
        writetable(temp.table, fileName, Sheet=jointsOfInterst{ang})
    end

end% end number of timepoints

