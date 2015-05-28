%% Lesion Volume for Local Directory
%
%   This script processes all vtu files in the chosen directory using
%   IsoAreaCompute to compute lesion volumes. A data structure is generated
%   containing the source file name and the computed lesion volume.
%
%   The main complication is the inclusion of the changing cell death
%   threshold. This is dealth with by loading experiment-#.mat

clear all
close all

% fp = ['../../RFA-sensitivity-cool/20150217-170854-experiment-1/run'];
% fp = ['../../RFA-sensitivity-cool/20150218-172328-experiment-2/run'];
% fp = ['../../RFA-sensitivity-cool/20150219-111815-experiment-3/run'];
% fp = ['../../RFA-sensitivity-cool/20150219-210431-experiment-4/run'];
% fp = ['../../RFA-sensitivity-cool/20150222-155751-experiment-5/run'];
fp = ['../../RFA-sensitivity-cool/20150223-120141-experiment-6/run'];

% fp_thresh = ['../../RFA-sensitivity-cool/experiment-1.mat'];
% fp_thresh = ['../../RFA-sensitivity-cool/experiment-2.mat'];
% fp_thresh = ['../../RFA-sensitivity-cool/experiment-3.mat'];
% fp_thresh = ['../../RFA-sensitivity-cool/experiment-4.mat'];
% fp_thresh = ['../../RFA-sensitivity-cool/experiment-5.mat'];
fp_thresh = ['../../RFA-sensitivity-cool/experiment-6.mat'];

fi = 'redundant'; % set field of interest  (this is currently redundant)
% thresh = 323; % set threshold for isoarea
% stem = 'enthalpy000'; % choose file stem -------------------------------------------------
thresh = 0.8; % set threshold for isoarea, this is also redundant now
stem = 'cell-death000'; % choose file stem -------------------------------------------------

% obtain file list from directory
file_list = dir([fp '*']);
file_list = file_list([file_list.isdir]); % get rid of anything that isnt a directory

% choose file numbers
j=[0:597];
% j=[0:100];

% load experiment
eval(['load ' fp_thresh])

% loop through time steps
for n=1:length(j)
    
    % create file stems
    s = length(file_list);
    for i=1:s
        file_list(i).name = [num2str(i) '/' stem sprintf('%03u',j(n)) '.vtu'];
    end

    % allocate struct
    results = struct('name',{file_list.name},'area',0.,'row',0.,'thresh',0.);


    % loop through all vtu files computing area 
    for i=1:s

        display(file_list(i).name)
        
        % set threshold
        results(i).row = sscanf(results(i).name,...
            ['%f/' stem '*']);
        results(i).thresh = thresh; % choose whether to adjust thresh or not ----------------------------------------------
%         results(i).thresh = (0.9-0.7)*A(results(i).row,12)+0.7;
        
        % collect responses
        [status, cmdout] = unix(['./IsoAreaCompute.py '...
            fp file_list(i).name ' ' num2str(results(i).thresh) ' ' fi])
        load lesion_area.mat;
        results(i).area = lesion_area;
        lesion_area

    end
    
    [a b] = sort([results.row]);
    results = results(b);
    
    % save
    full_results(n,:)=results;
    
    clear results

end


eval(['save ' fp 'lesion_area.mat full_results j'])