%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function test_rock_filename = text_to_matlab
%
% This function is designed to be called by lasergui.m. It shall prompt the
% user for a directory containing text files that they wish to test with
% their PLS model. The function will look to see if the folder contains a
% .mat file containing the data already converted from .txt. If this .mat
% file exists, the function returns the file. Otherwise it will convert
% each txt file to a vector and concatenate them to the test_rock_data
% matrix (one txt file's observations per row), then save and return this
% matrix.
%
% This version of the code is meant to work only with test data 12288 in
% length.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function test_rock_filename = text_to_matlab(dir_in,mode)

% Save original directory and prompt user for directory containing test
% data
orig_dir = pwd;
test_rock_filename = [];

switch mode
    case 1
        singlefolder = dir_in;
        setcount = 1;
        setfoldername = 'NA';
    case 2
        setfolder = dir_in;
        [basedir, setfoldername] = fileparts(setfolder);
        cd(setfolder)
        setlist = dir;
        setlist = extractfield(setlist,'name');
        setlist = setlist(3:numel(setlist));
        setcount = numel(setlist); 
end

w = waitbar(0,'Converting testing data to .mat format...');
date_n_time = datestr(now);
colon_find = strfind(date_n_time,':');
date_n_time(colon_find) = '_';
waitbar(1/5,w)
for t = 1:setcount
    waitbar((1/5+3/5*(t/setcount)),w)
    % Determine folder if in Testing Set mode
    switch mode
        case 1
            folder = singlefolder;
        case 2
            folder = [basedir,'\',setfoldername,'\',setlist{t}];
    end
    % Use user-specified folder to generate filename of .mat file
    cd(folder)
    cd ..
    dir_length = length(pwd) + 2;
    cd(folder)
    rock_type = folder(dir_length:length(folder));
    test_rock_filename{t} = strcat(rock_type,'.mat');

    % Get the filenames for each .txt file in a list and the number of files.
    dirListing = dir(fullfile(folder,'*.txt'));
    nfiles = length(dirListing);

    % Initialize test_rock_data matrix.
    test_rock_data = [];

    % Place data from each text file on its own row in test_rock_data
    for i=1:nfiles
        fileName = fullfile(folder,dirListing(i).name);
        [~,name,~] = fileparts(fileName);
        test_rock_data(i,:) = dlmread(fileName,'\t','B11..B12298');
    end
    
    waitbar(1,w)
    temp_dir = pwd;
    save_dir = check_create_dir(['Results\Testing Data - Conversion to mat\',date_n_time],mfilename('fullpath'),2);
    cd(save_dir)
    save(test_rock_filename{t}, 'test_rock_data');
    cd(temp_dir)
    disp(['Converted test data for ', rock_type, ' saved to ', save_dir])
    cd ..

end
cd(orig_dir)
delete(w)
