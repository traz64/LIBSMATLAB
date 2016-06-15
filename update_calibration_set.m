%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function update_calibration_set(calibration_dir, XRF_spreadsheet, save_dir)
%
% This function updates the calibration set to include all of the XRF data
% in the spreadsheet specified by the input to this function. It will
% create and save a struct to the specified directory that contains
% spectrometer data and corresponding XRF data along with class labels
% specifying the rock type for each observation (rows) so that we can
% ensure that the XRF data directly corresponds to the spectrometer data.
%
% Before running this function, make sure that the spectrometer data for
% each type of rock has been converted to a .txt file and stored in
% 'calibration_dir'. The spectrometer data files should be named with the
% following convention: MM-DD-YYYY Rock Name Sam#Loc#.txt
% For example, if the spectrometer data was taken on July 22, 2015 and the
% rock type was Atkinson Quartzite, using the first sample and the third
% location, the filename would be:
%                       07-22-2015 Atkinson Quartzite Sam1Loc3.txt
%
% Inputs:
%   calibration_dir - directory that contains the spectrometer data to be
%   used for calibration of the PLS system.
%
%   XRF_spreadsheet - filename of spreadsheet containing XRF data including
%   direct path to file location. File location must be included in the
%   MATLAB path.
%
%   save_dir - directory to which the calibration set will be saved.
%
% Outputs:
%   Calibration set - A structure containing the spectrometer data and
%   corresponding XRF data, to be saved to the specified directory.
%   Filename will be 'DD-MMM-YYYY Calibration Set.mat' (e.g. 22-Jul-2015
%   Calibration Set.mat'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function update_calibration_set(calibration_dir, XRF_spreadsheet, save_dir)

% Save original directory to switch back to it after calibration set has
% been saved.
original_dir = pwd;

% Determine number of rock types to be used in calibration.
dirListing = dir(fullfile(calibration_dir,'*.txt'));
nfiles = (length(dirListing))/50;

% Read XRF spreadsheet to obtain list of rock types to be used in
% calibration set.
cells = ['A1:Z',num2str(nfiles+1)];
[~, ~, XRF_data] = xlsread(XRF_spreadsheet,1,cells);
XRF_labels = XRF_data(1:nfiles,1);

% Initialize X matrix (will contain spectrometer data, 1 observation per
% row)
X = zeros(length(dirListing),12288);

% Place data from each text file on its own row in test_rock_data. Also
% populate X_labels matrix with rock type of each observation.
for i=1:nfiles*50
    fileName = fullfile(calibration_dir,dirListing(i).name);
    rockType = fileName(length(calibration_dir)+12:length(fileName));
    rockType = fliplr(rockType);
    rockType(1:13) = [];
    rockType = fliplr(rockType);
    X_labels{i} = strtrim(rockType);
    X(i,:) = dlmread(fileName,'\t','B11..B12298');
end
X_labels = X_labels';

% Initialize Y matrix
Y = [];

% Create Y matrix containing XRF data corresponding to the X_labels matrix.
for i = 1:size(X_labels,1)
    label_index = find(strcmp(X_labels{i},XRF_labels));
    if isempty(label_index)
        disp(['error on rock ',X_labels{i}])
        disp('not in list')
        disp(XRF_labels)
   
    else disp('rock ok')
    end
    C(i,1)=XRF_data{label_index,2};
    for j = 3:size(XRF_data,2)
        Y(i,j-2) = XRF_data{label_index,j};
    end
end

% Place X, Y, and X_labels in a struct and save in specified directory.
rock_data = struct('X',X,'Y',Y,'C',C);
if sum(isnan(Y))
    disp('not gonna work')
end
cd(save_dir)
save([date,' Calibration Set.mat'],'rock_data')

% Change back to original directory and output completion message.
cd(original_dir)
disp('Finished updating calibration set.');
disp([date,' Calibration Set.mat saved to ', save_dir])
disp('')
