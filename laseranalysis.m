%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function laseranalysis(rock_data, PLS, ncomp_in, mode_custom_thresh)
%
% This function handles the training and testing of a PLS model to predict
% the chemical components of an observed data set. In training mode, this
% function will take the observed X data and known Y data as well as the
% number of PLS components to use and train a PLS model, saving the beta
% matrix for use in testing. In testing mode, the function loads the saved
% beta matrix and creates a predicted chemical composition for the observed
% data set.
%
% Inputs:
%   rock_data - struct containing X and Y data with same number of observations
%   and corresponding rows aligned (training) or a matrix containing all of
%   the X data (testing)
%       (X - observed data set. In training this contains the spectrometer data
%       for all aggregate stones. In testing, this contains the spectrometer
%       data for just one aggregate stone.
%
%       Y - known data set. In training, this contains the known composition of
%       all aggregate stones. For testing, enter 0 for this input.)
%
%   PLS - struct containing PLS generated in training mode. Only used in
%   testing mode. For training mode, this input will be null.
%
%   ncomp - number of PLS components to be used in training the PLS model.
%   In training mode enter the number of PLS components to use or enter a
%   number <= 0 to automatically use the maximum PLS components. In testing
%   mode, enter 0 for this input.
%
%   mode - 'train' or 'test' (obtained from GUI)
%
% Outputs:
%   Training mode - Saves the PLS model in a struct; contains beta matrix
%   necessary for testing.
%   Testing mode - Saves the predicted chemical composition for the input
%   observed spectrometer data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function laseranalysis(rock_data, PLS, ncomp_in, mode, thresholds)
%% TRAINING MODE %%
if(strcmp(mode,'train'))
    w = waitbar(0,'Generating PLS Model...');
    current_dir = check_create_dir('void',mfilename('fullpath'),0);
    
    % LOAD DATA IF FUNCTION IS CALLED FROM GUI
    if(ischar(rock_data))
        disp('Loading Rock Data.')
        rock_data = load(rock_data);
        rock_data = rock_data.rock_data;
        X = rock_data.X;
        
        disp('Loading Y Data.')
        Y = rock_data.Y;
        classer = rock_data.C;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % PREPROCESSING STAGE %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NOTE %
    %%%%%%%%
    % This stage can eventually contain all of the different types of
    % preprocessing techniques we want to test. This includes the split
    % training technique, Y-scaling, and normalizing the data to total
    % light emission.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    waitbar(1/10,w);
    % Remove light intensity values less than 0.
    [m,n]=size(X);

    disp('Removing negative light intensity values.')
    for i=1:m
        for j=1:n
            if X(i,j)<0
                X(i,j)=0;
            end
        end
    end
    waitbar(3/10,w);
    % Normalizing Data to Total Light Emission
    disp('Normalizing spectra to total light emissions.')
    
    % Initializing total light intensity.
        total_light_int=zeros(m,1);
    
    for i=1:m
        for j=1:n
            total_light_int(i) = total_light_int(i) + X(i,j);
        end
    end
    
    % Determine Xnorm
    Xnorm = zeros(m,n);
    for i=1:m
        for j=1:n
            Xnorm(i,j) = X(i,j)/total_light_int(i);
        end
    end
    waitbar(4/10,w);
    % Perform Y Scaling
    numCol= size(Y,2);
    
    % Initialize minmax matrix. Currently 2x24, but it can be of varying
    % size.
    minmax_base = zeros(2, numCol);
    for i = 1:numCol
        
        maxVal = max(Y(:,i));
        minVal = min(Y(:,i));
        minmax_base(2,i) = minVal;
        minmax_base(1,i) = maxVal;
        val_range = maxVal-minVal;
        if(val_range == 0)
            Yscaled(:,i) = 0;
        else
            Yscaled(:,i) = (Y(:,i) - minVal)/val_range;
        end
        
    end
    waitbar(5/10,w);
    % NOTE: Do Y scaling for Carb and Nonc models after splitting into
    % separate matrices
    
    % Determine pls components if set to 'Auto'
    % THIS NEEDS TO BE DONE SEPARATELY FOR EACH OF THE SPLIT TRAINING
    % MODELS
    
    ncomp_base_in = ncomp_in.Base;
    if ncomp_base_in<=0;
        disp('Determining optimal number of PLS components for Base model.')
        ncomp = 25;
        
        % Find PLS regression for starting number of ncomp
        [~,~,~,~,~,PCTVAR] = plsregress(Xnorm,Yscaled,ncomp);
        pctvar_count = 0;
        PCTVAR = 100*PCTVAR;
        
        % Set ncomp to max number that explains >1% of variation
        for i = 1:ncomp
            if (PCTVAR(2,i) >= 1)
                pctvar_count = pctvar_count + 1;
            else
                break;
            end
        end
        ncomp = i - 1;
        disp(['Number of Base PLS components automatically set to ', num2str(ncomp),'.']);
    else
        ncomp = ncomp_base_in;
    end
    waitbar(7/10,w);
    % PLS REGRESSION
    disp('Generating Broad Base PLS Model.')
    [~,~,~,~,betamat,PCTVAR] = plsregress(Xnorm,Yscaled,ncomp);
    
    % Creates a structure with all of the PLS Model variables to be saved as a
    % .mat file.
    PLS_Model_base = struct('Date', date, 'PercentVariation', PCTVAR, 'Beta', betamat, 'MinMax', minmax_base);
    PLS_Model_base.NComp = ncomp;
    
    % Save PLS Model in results directory and return to pwd.
    dir = pwd;
    [~] = check_create_dir('PLS models',mfilename('fullpath'),1);
    cd([current_dir,'\PLS models'])
    cd(dir);
    
    % SPLIT TRAINING
    Ycarb = [];
    Xcarb = [];
    Ynonc = [];
    Xnonc = [];    
    Ytrap = [];
    Xtrap = [];
    for i = 1:size(Y,1)
        switch classer(i)
            
            case 'c'  
                Ycarb = cat(1,Ycarb,Y(i,:));
                Xcarb = cat(1,Xcarb,Xnorm(i,:));
            case 'n'
                Ynonc = cat(1,Ynonc,Y(i,:));
                Xnonc = cat(1,Xnonc,Xnorm(i,:));
            case 't'
                Ytrap = cat(1,Ytrap,Y(i,:));
                Xtrap = cat(1,Xtrap,Xnorm(i,:));
            otherwise
                error('Error in XRF Sheet used in Calibration. Unknown classification')
        end
    end
    
    % Calculate ncomp for carbonate rocks
    ncomp_carb_in = ncomp_in.Carbonate;
    if ncomp_carb_in<=0;
        disp('Determining optimal number of PLS components for Carbonate model.')
        ncomp_carb = 25;
        
        % Find PLS regression for maximum number of ncomp
        [~,~,~,~,~,PCTVAR] = plsregress(Xcarb,Ycarb,ncomp_carb);
        pctvar_count = 0;
        PCTVAR = 100*PCTVAR;
        % Set ncomp to max number that explains >1% of variation
        % Stop searching if two components in a row are found to explain
        % less than 1%.
        for i = 1:ncomp_carb
            if (PCTVAR(2,i) >= 1)
                pctvar_count = pctvar_count + 1;
            else
                break;
            end
        end
        ncomp_carb = i - 1;
        disp(['Number of Carbonate PLS Model components automatically set to ', num2str(ncomp_carb),'.']);
    else
        ncomp_carb = ncomp_carb_in;
    end
    
    % Perform Y scaling for Carbonate rocks
    minmax_carb = zeros(2, numCol);
    for i = 1:numCol
        
        maxVal = max(Ycarb(:,i));
        minVal = min(Ycarb(:,i));
        minmax_carb(2,i) = minVal;
        minmax_carb(1,i) = maxVal;
        val_range = maxVal-minVal;
        if(val_range == 0)
            Yscaled_carb(:,i) = 0;
        else
            Yscaled_carb(:,i) = (Ycarb(:,i) - minVal)/val_range;
        end
        
    end
    %-------------------------------------------------------------
    % Calculate ncomp for noncarbonate rocks
    ncomp_nonc_in = ncomp_in.Carbonate;
    if ncomp_nonc_in<=0;
        disp('Determining optimal number of PLS components for Non-Carbonate model.')
        ncomp_nonc = 25;
        
        % Find PLS regression for maximum number of ncomp
        [~,~,~,~,~,PCTVAR] = plsregress(Xnonc,Ynonc,ncomp_nonc);
        pctvar_count = 0;
        PCTVAR = 100*PCTVAR;
        % Set ncomp to max number that explains >1% of variation
        % Stop searching if two components in a row are found to explain
        % less than 1%.
        for i = 1:ncomp_nonc
            if (PCTVAR(2,i) >= 1)
                pctvar_count = pctvar_count + 1;
            else
                break;
            end
        end
        ncomp_nonc = i - 1;
        disp(['Number of Non-Carbonate PLS Model components automatically set to ', num2str(ncomp_nonc),'.']);
    else
        ncomp_nonc = ncomp_nonc_in;
    end
    
    % Perform Y scaling for NonCarbonate rocks
    minmax_nonc = zeros(2, numCol);
    for i = 1:numCol
        
        maxVal = max(Ynonc(:,i));
        minVal = min(Ynonc(:,i));
        minmax_nonc(2,i) = minVal;
        minmax_nonc(1,i) = maxVal;
        val_range = maxVal-minVal;
        if(val_range == 0)
            Yscaled_nonc(:,i) = 0;
        else
            Yscaled_nonc(:,i) = (Ynonc(:,i) - minVal)/val_range;
        end
        
    end
    %------------------------------------------------------------
    ncomp_trap_in = ncomp_in.Trap;
    if ncomp_trap_in<=0;
        disp('Determining optimal number of PLS components for Traprock model.')
        ncomp_trap = 25;
        
        % Find PLS regression for maximum number of ncomp
        [~,~,~,~,~,PCTVAR] = plsregress(Xcarb,Ycarb,ncomp_trap);
        pctvar_count = 0;
        PCTVAR = 100*PCTVAR;
        % Set ncomp to max number that explains >1% of variation
        % Stop searching if two components in a row are found to explain
        % less than 1%.
        for i = 1:ncomp_trap
            if (PCTVAR(2,i) >= 1)
                pctvar_count = pctvar_count + 1;
            else
                break;
            end
        end
        ncomp_trap = i - 1;
        disp(['Number of Traprock PLS Model components automatically set to ', num2str(ncomp_trap),'.']);
    else
        ncomp_trap = ncomp_trap_in;
    end
    
    % Perform Y scaling for Traprocks
    minmax_trap = zeros(2, numCol);
    for i = 1:numCol
        
        maxVal = max(Ytrap(:,i));
        minVal = min(Ytrap(:,i));
        minmax_trap(2,i) = minVal;
        minmax_trap(1,i) = maxVal;
        val_range = maxVal-minVal;
        if(val_range == 0)
            Yscaled_trap(:,i) = 0;
        else
            Yscaled_trap(:,i) = (Ytrap(:,i) - minVal)/val_range;
        end
        
    end

    
    waitbar(9/10,w);
    
    % PLS REGRESSION Trap
    disp('Generating Trap Rock PLS Model.')
    [~,~,~,~,betamat,PCTVAR] = plsregress(Xtrap,Yscaled_trap,ncomp_trap);
    
    % Creates a structure with all of the PLS Model variables to be saved as a
    % .mat file.
    PLS_Model_trap = struct('Date', date, 'PercentVariation', PCTVAR, 'Beta', betamat, 'MinMax', minmax_trap);
    PLS_Model_trap.NComp = ncomp_trap;
    
    % Save PLS Model in results directory and return to pwd.
    dir = pwd;
    cd([current_dir,'\PLS models'])
    cd(dir);  
    
    % PLS REGRESSION CARB
    disp('Generating Carbonate Rock PLS Model.')
    [~,~,~,~,betamat,PCTVAR] = plsregress(Xcarb,Yscaled_carb,ncomp_carb);
    
    % Creates a structure with all of the PLS Model variables to be saved as a
    % .mat file.
    PLS_Model_carb = struct('Date', date, 'PercentVariation', PCTVAR, 'Beta', betamat, 'MinMax', minmax_carb);
    PLS_Model_carb.NComp = ncomp_carb;
    
    % PLS REGRESSION NONC
    disp('Generating Non-Carbonate Rock PLS Model.')
    [~,~,~,~,betamat,PCTVAR] = plsregress(Xnonc,Yscaled_nonc,ncomp_nonc);
    
    % Creates a structure with all of the PLS Model variables to be saved as a
    % .mat file.
    PLS_Model_nonc = struct('Date', date, 'PercentVariation', PCTVAR, 'Beta', betamat, 'MinMax', minmax_nonc);
    PLS_Model_nonc.NComp = ncomp_nonc;
    
    % Save PLS Model in results directory and return to pwd.
    dir = pwd;
    date_n_time = datestr(now);
    colon_find = strfind(date_n_time,':');
    date_n_time(colon_find) = '_';
    cd([current_dir,'\PLS Models'])
    
    waitbar(10/10,w);
    
    % Create PLS structure that contains PLS models for base, carbonate,
    % and noncarbonate all together and save. This is the file the user
    % should load when testing the system as it contains all of the
    % necessary PLS models.
    PLS_Model_All = struct('Base', PLS_Model_base, 'Carbonate', PLS_Model_carb, 'NonCarbonate', PLS_Model_nonc, 'Trap', PLS_Model_trap);
    save(['PLS-Model-All-', date_n_time, '.mat'], 'PLS_Model_All');
    cd(dir);
    disp(['PLS Model saved to ', current_dir, '\PLS Models'])
    disp('Model calibration complete.')
    delete(w)
end
%% TESTING MODE %%
close all;
if(strcmp(mode,'test')) || (strcmp(mode, 'testset'))
    w = waitbar(0,'Processing testing data...');
    setdata = rock_data;
    sampnum = numel(setdata);
    date_n_time = datestr(now);
    colon_find = strfind(date_n_time,':');
    date_n_time(colon_find) = '_';
    current_dir = check_create_dir('void',mfilename('fullpath'),0);
    save_dir = check_create_dir(['Results\Testing Data - Analysis\',date_n_time],mfilename('fullpath'),2);
    means=cell(sampnum+1,25);
    chem = {'SiO2' 'Al2O3' 'Fe2O3' 'CaO' 'MgO' 'Na2O' 'P2O5' 'TiO2' 'K2O' 'MnO' 'BaO' 'SO3' 'SrO' 'CuO' 'ZrO2' 'ZnO' 'Y2O3' 'Rb2O' 'Ga2O3' 'Cl' 'Cr2O3' 'NiO' 'CeO2' 'Nb2O5'};
    means(1,2:25)=chem;
    stddev=means;
    for c = 1:sampnum
        waitbar((1/15+14/15*(c/sampnum)),w)
        rock_data = setdata{c};
        % LOAD BETA MATRIX    
        dir = pwd;
        disp([setdata{c},':'])
        disp('Loading PLS model.')
        if(ischar(PLS))
            load(PLS);
        end
        betamat = PLS_Model_All.Base.Beta;

        % Obtain full path to rock_data file
        rock_data_location = which(rock_data);
        rock_data_location = rock_data_location(1:length(rock_data_location)-length(rock_data));

        % Obtain name of rock from rock_data filename.
        rock_type = rock_data(1:length(rock_data)-4);

        if(ischar(rock_data))
            disp('Loading X Data.')
            rock_data = load(rock_data);
            X = rock_data.test_rock_data;
        end

        % Preprocessing stage
        % Should use the same techniques as in training.


        [m,n]=size(X);
        % Center Clipping
        for i=1:m
            for j=1:n
                if X(i,j)<0
                    X(i,j)=0;
                end
            end
        end

        % Normalizing Data to Total Light Emission
        disp('Normalizing spectra to total light emission.')

        % Initializing total light intensity.

            total_light_int = zeros(1,m);

        for i=1:m
            for j=1:n
                total_light_int(i) = total_light_int(i) + X(i,j);
            end
        end

        % Initialize Xnorm
        Xnorm = zeros(m,n);
        for i=1:m
            for j=1:n
                Xnorm(i,j) = X(i,j)/total_light_int(i);
            end
        end


        % MAKE INITIAL PREDICTION
        Ypredicted = [ones(size(Xnorm,1),1) Xnorm]*betamat;

        % Obtain beta and minmax from carbonate or non-carbonate models based
        % on carbonate threshold. Beta is used for prediction, minmax is used
        % for reverse Y scaling.

        % Edit: Reverse Y scaling has been changed back to just using the Base
        % min_max values rather than split based on Carbonate content.
        min_max = PLS_Model_All.Base.MinMax;

        numCol = size(Ypredicted,2);
        for i = 1:numCol
            val_range = min_max(1,i)-min_max(2,i);
            Ypredicted(:,i) = (Ypredicted(:,i)*val_range) + min_max(2,i);
        end

        % Use prediction matrix to determine whether the rock is carbonate,
        % non-carbonate, or trap, and then make another prediction using the
        % corresponding PLS Model.
        
        % Ratio for classification
        threshratio = mean((Ypredicted(:,1)./abs(Ypredicted(:,4))).^2.*abs(Ypredicted(:,3)));
        
        if (thresholds == (-1 -1))
            carbthresh=150;
            noncarbthresh=1150;
        else
            carbthresh=thresholds(1);
            noncarbthresh=thresholds(2);
        end
           
        if(threshratio <= carbthresh)
             betamat = PLS_Model_All.Carbonate.Beta;
             min_max = PLS_Model_All.Carbonate.MinMax;
             disp('Classified as Carbonate')
        elseif (threshratio<=noncarbthresh)
            betamat = PLS_Model_All.Trap.Beta;
            min_max = PLS_Model_All.Trap.MinMax;
            disp('Classified as Trap')
        else
            betamat = PLS_Model_All.NonCarbonate.Beta;
            min_max = PLS_Model_All.NonCarbonate.MinMax;
            disp('Classified as Non-Carbonate')
        end

        % Make new prediction based on split training decision.
        Ypredicted = [ones(size(Xnorm, 1), 1) Xnorm] * betamat;

        % Perform reverse Y scaling on predicted matrix.
        numCol = size(Ypredicted,2);
        for i = 1:numCol

            val_range = min_max(1,i)-min_max(2,i);
            %if(range == 0)
            %   Ypredicted(:,i) = 0;
            %else
            Ypredicted(:,i) = (Ypredicted(:,i)*val_range) + min_max(2,i);
            %end

        end

        % MEAN AND STD. DEV. CALCULATIONS
        % Initialize mean and std. dev. variables
        Ymean = zeros(1,size(Ypredicted,2));
        Ystd = zeros(1,size(Ypredicted,2));
        Ymode = zeros(1,size(Ypredicted,2));
        Ymedian = zeros(1,size(Ypredicted,2));

        % Calculate mean and std. dev. for each column.
        for i = 1:size(Ypredicted,2)
            Ymean(i) = mean(Ypredicted(:,i));
            Ystd(i) = std(Ypredicted(:,i));
            Ymedian(i) = median(Ypredicted(:,i));
        end

        % DISPLAY RESULTS IN FORMATTED TABLE
        figure('Name',rock_type,'NumberTitle','off')
        Ymstd = cat(1,Ymean,Ystd,Ymedian);
        title(['Results - ' rock_type])
        t = uitable('Data', Ypredicted, 'ColumnName', chem);
        set(t,'Position',[0 80 560 315])
        h = uitable('Data', Ymstd, 'ColumnName', chem, 'RowName', {'Mean','Std.','Median'});
        set(h,'Position',[0 0 560 80])
        %data extraction for excel
            rockname=setdata{c};
            means(c+1,1)= {rockname(1:end-4)};
            stddev(c+1,1)= {rockname(1:end-4)};
            
            means(c+1,2:25)= num2cell(Ymean);
            stddev(c+1,2:25)= num2cell(Ystd);
            
        % COMPILE RESULTS STRUCTURE
        Y_Results = struct('Y_Predicted', Ypredicted, 'Y_Mean', Ymean, 'Y_Std', Ystd, 'Y_Median', Ymedian);

        % SAVE PREDICTION MATRIX
        cd(save_dir)
        save([rock_type,' Analysis Results.mat'], 'Y_Results');
        cd(dir);
        disp(['Prediction saved to ', save_dir])
        fprintf('\n')
    end
    if exist('actxserver','file')
        if(strcmp(mode,'test'))
            ExcelName = 'Single Test';
        elseif(strcmp(mode,'testset'))
            ExcelName = 'Testing Set';
        end
        cd(save_dir)
        warning('off','MATLAB:xlswrite:AddSheet')
        xlswrite(['{Mean & Std}' ExcelName ' Results Summary.xls'],means,'Mean Values')
        xlswrite(['{Mean & Std}' ExcelName ' Results Summary.xls'],stddev,'Standard Deviations')
        disp(['Prediction report spreadsheet saved to ', save_dir])
        fprintf('\n')
        objExcel = actxserver('Excel.Application');
        objExcel.Workbooks.Open(fullfile(cd,['{Mean & Std}' ExcelName ' Results Summary.xls'])); % Full path is necessary!
        try
          objExcel.ActiveWorkbook.Worksheets.Item(1).Delete;
          invoke(objExcel.ActiveWorkbook.Worksheets.Item(1), 'Activate');
        catch
        end
        % Save, close and clean up.
        objExcel.ActiveWorkbook.Save;
        objExcel.ActiveWorkbook.Close;
        objExcel.Quit;
        objExcel.delete;
        msgbox(['An Excel spreadsheet containing a summary of this test has been saved in:' save_dir])
    else
        warning('ActiveX process could not be created, excel summary not saved ')
    end
    cd(dir);
    delete(w)
end