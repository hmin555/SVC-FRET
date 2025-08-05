%%% ====================== MAIN SCRIPT ====================== %%%
% SVC-FRET Calibration System
% Author: [Min Hu]
% Date: 2025-08-04
% Description: Performs FRET system calibration using SVC-FRET method with
%              three standard plasmids of known FRET efficiencies

% --- Environment Initialization ---
clear;          % Clear workspace variables
close all;       % Close all open figures
clc;             % Clear command window

% --- FRET Image Path Configuration ---
% Define paths to plasmid image directories (C4Y, C10Y, C40Y, C80Y)
C4Y_path = "H:\SVC-FRET代码及数据\01标准质粒测试 60X OD2\Test\C4Y";
C10Y_path = "H:\SVC-FRET代码及数据\01标准质粒测试 60X OD2\Test\C10Y";
C40Y_path = "H:\SVC-FRET代码及数据\01标准质粒测试 60X OD2\Test\C40Y";
C80Y_path = "H:\SVC-FRET代码及数据\01标准质粒测试 60X OD2\Test\C80Y";
paths = [C4Y_path; C10Y_path; C40Y_path; C80Y_path];  % Path collection

% --- Data Structures Initialization ---
fittedVector = zeros(3);      % Stores fitted 3D vectors
I_DD = zeros(3, 1);           % Donor channel intensities
I_DA = zeros(3, 1);           % FRET channel intensities
I_AA = zeros(3, 1);           % Acceptor channel intensities

% Known FRET efficiencies for standard plasmids
standPlasmidEDs = [0.299, 0.223, 0.158];  % Donor efficiencies
standPlasmidEAs = [0.299, 0.223, 0.158];  % Acceptor efficiencies

% --- Processing Pipeline ---
tic;  % Start performance timer

%%% --- STEP 1: Load multi-field three-channel data ---
for i = 1:3
    field_data = process_fret_data(char(paths(i))); 
    
    % Store data by plasmid type
    switch i
        case 1
            allPointsC4Y = struct2cell(field_data);  
        case 2
            allPointsC10Y = struct2cell(field_data);       
        case 3
            allPointsC40Y = struct2cell(field_data);       
    end
end

%%% --- STEP 2: Random field selection and vector fitting ---
num_trials = 1000;     % Number of random trials
seed = 0;           % Random seed for reproducibility
%results = struct(); % Results storage

for trial = 1:num_trials
    %rng('shuffle');
    rng(trial + seed);  % Set reproducible random seed[1](@ref)
    
    % Randomly select one field per plasmid
    idx_C4Y = randi(length(allPointsC4Y));
    idx_C10Y = randi(length(allPointsC10Y));
    idx_C40Y = randi(length(allPointsC40Y));
    
    % Extract selected data
    selected_C4Y = allPointsC4Y{idx_C4Y};
    selected_C10Y = allPointsC10Y{idx_C10Y};
    selected_C40Y = allPointsC40Y{idx_C40Y};

    % Fit 3D vectors using PCA
    fittedVector(1, :) = fit_fret_vectors(selected_C4Y);
    fittedVector(2, :) = fit_fret_vectors(selected_C10Y);
    fittedVector(3, :) = fit_fret_vectors(selected_C40Y);
    
    [a,b]=evaluateCoplanarity(fittedVector(1, :),fittedVector(2, :),fittedVector(3, :),1e-3);
    %if a==1
    % Extract channel intensities
    for i = 1:3
        I_DD(i) = fittedVector(i, 1);
        I_DA(i) = fittedVector(i, 2);
        I_AA(i) = fittedVector(i, 3);
    end

    % Calculate system parameters
    [~, sysParams] = calculate_system_parameters(...
        I_DD, I_DA, I_AA, standPlasmidEDs, standPlasmidEAs);
    
    % Store additional metadata (filed from plasmid)
    sysParams.C4Y = allPointsC4Y{idx_C4Y}(1, 6);
    sysParams.C10Y = allPointsC10Y{idx_C10Y}(1, 6);
    sysParams.C40Y = allPointsC40Y{idx_C40Y}(1, 6);

    % Aggregate results
    results(trial) = sysParams;
    %end
end

% --- Results Processing ---
dataTable = struct2table(results);  % Convert to table

% Statistical analysis (mean ± SD)
keyParams = {'G', 'a', 'd', 'k', 'eA_eD'};  % Critical parameters
statsReport = table();
valid_idx = true(height(dataTable),1);
for param = keyParams
    paramName = param{1};
    mu = mean(dataTable.(paramName));
    sigma = std(dataTable.(paramName));
    statsReport.(paramName) = [mu, sigma];
    valid_idx = valid_idx & (dataTable.(param{1}) >= mu-1*sigma) & (dataTable.(param{1}) <= mu+1*sigma);
end
filtered_data = dataTable(valid_idx, :);
allMeans = varfun(@mean, filtered_data);
disp('System Parameters Summary:');
disp(statsReport);

toc;  % Display elapsed time
%%% ==================== END MAIN SCRIPT ===================== %%%

%%% ====================== FUNCTIONS ======================== %%%

function fretData = read_fret_images(mainDir)
% READ_FRET_IMAGES Loads FRET image data from subdirectories
% Input:  mainDir - Main directory containing FRET_XXXX subfolders
% Output: fretData - Structured array with image data and metadata
%
% Features:
% - Automatically detects FRET_XXXX subfolders
% - Handles multi-page TIFF stacks
% - Validates folder naming conventions[3](@ref)

% Validate input directory
if ~isfolder(mainDir)
    error('Directory not found: %s', mainDir);
end

% Find FRET_XXXX subdirectories
subFolders = dir(fullfile(mainDir, 'FRET_*'));
subFolders = subFolders([subFolders.isdir]);

if isempty(subFolders)
    error('No valid FRET_* subfolders found');
end

% Initialize output structure
fretData = struct('folderName', {}, 'fourDigits', {}, ...
                  'AA', {}, 'DA', {}, 'DD', {});

% Process each subfolder
for i = 1:length(subFolders)
    folderName = subFolders(i).name;
    fourDigits = folderName(6:9);  % Extract numerical identifier
    
    % Validate folder name format
    if ~all(isstrprop(fourDigits, 'digit'))
        warning('Skipping invalid folder: %s', folderName);
        continue;
    end
    
    % Build full TIFF paths
    subPath = fullfile(mainDir, folderName);
    fileNames = {
        fullfile(subPath, [fourDigits 'AA.tif']), ...
        fullfile(subPath, [fourDigits 'DA.tif']), ...
        fullfile(subPath, [fourDigits 'DD.tif'])
    };
    
    % Read and store image data
    imageData = cell(1, 3);
    for j = 1:3
        if ~isfile(fileNames{j})
            error('Image not found: %s', fileNames{j});
        end
        
        info = imfinfo(fileNames{j});
        numPages = numel(info);
        firstPage = imread(fileNames{j}, 1);
        
        % Handle multi-page TIFFs
        if numPages > 1
            imgStack = zeros([size(firstPage), numPages], 'like', firstPage);
            for k = 1:numPages
                imgStack(:, :, k) = imread(fileNames{j}, k);
            end
            imageData{j} = imgStack;
        else
            imageData{j} = imread(fileNames{j});
        end
    end
    
    % Populate output structure
    fretData(i).folderName = folderName;
    fretData(i).fourDigits = fourDigits;
    fretData(i).AA = imageData{1};  % Acceptor excitation
    fretData(i).DA = imageData{2};  % FRET channel
    fretData(i).DD = imageData{3};  % Donor excitation
end

fprintf('Processed %d subfolders\n', length(fretData));
end

function allFOV = process_fret_data(imagesDir)
% PROCESS_FRET_DATA Processes FRET images to extract valid data points
% Input:  imagesDir - Directory containing FRET images
% Output: allFOV - Structure with processed data per field-of-view
%
% Processing steps:
% 1. Background subtraction using histogram-based estimation
% 2. Foreground masking via intensity thresholding
% 3. Field flattening using circular mask
% 4. Data pooling for dimensionality reduction)

% Load raw image data
fretData = read_fret_images(imagesDir);
numFields = numel(fretData);

% Initialize output
imagePoints = [];

for i = 1:numFields
    fieldID = str2double(fretData(i).fourDigits);
    imgSize = size(fretData(i).AA);
    
    % Initialize image stacks
    rawStack = cat(3, fretData(i).AA, fretData(i).DA, fretData(i).DD);
    rawStack = double(rawStack);  % Convert to double precision
    fgStack = zeros(size(rawStack)); % Foreground stack
    bgStack = zeros(size(rawStack));  % Background stack
    
    % --- Background Estimation ---
    for j = 1:3
        channelData = rawStack(:, :, j);
        
        % Estimate background from intensity histogram
        validPixels = channelData(channelData > 0 & channelData < 200);
        [counts, bins] = histcounts(validPixels, 0:200);
        [~, maxIdx] = max(counts);
        bgValue = bins(maxIdx);
        
        % Create background image and subtract
        bgStack(:, :, j) = bgValue;
        fgStack(:, :, j) = channelData - bgValue;
    end
    % --- Foreground Masking ---
    maskStack = rawStack - 3 * bgStack;  % Create mask
    maskStack(maskStack < 0) = 0;        % Threshold low values
    maskStack(maskStack > 65535) = 0;    % Clip overflow
    maskStack(maskStack > max(bgStack(:)) * 5) = 0; % Remove bright outliers    

    % Combine masks across channels
    combinedMask = all(maskStack > 0, 3);
       
    % --- Field Flattening ---
    center = (imgSize + 1) / 2;
    [X, Y] = meshgrid(1:imgSize(2), 1:imgSize(1));
    distances = sqrt((X - center(2)).^2 + (Y - center(1)).^2);
    fieldMask = double(distances <= sqrt(imgSize(1)^2 * 0.85 / pi));
    
    % Apply masks and extract valid pixels
    validPixels = fgStack .* combinedMask.* fieldMask;

    % --- Data Pooling ---
    pooledSize = 2048;  % Reduced resolution for efficiency
    pooledStack = zeros(pooledSize, pooledSize, 3);
    for j = 1:3
        pooledStack(:, :, j) = imresize(validPixels(:, :, j), [pooledSize pooledSize]);
    end
    
    % Extract coordinates of valid pixels
    validMask = all(pooledStack > 0, 3);
    [row, col] = find(validMask);
    pixelX = pooledStack(sub2ind(size(pooledStack), row, col, 3*ones(size(row)))); % DD
    pixelZ = pooledStack(sub2ind(size(pooledStack), row, col, ones(size(row)))); % DA
    pixelY = pooledStack(sub2ind(size(pooledStack), row, col, 2*ones(size(row)))); % AA

    
    % Compose data matrix
    imagePoints = [imagePoints; col, row, pixelX, pixelY, pixelZ, repmat(fieldID, numel(row), 1)];
end

% Organize by field-of-view
fovIDs = unique(imagePoints(:, 6));
allFOV = struct();
for k = 1:length(fovIDs)
    fovName = ['FOV_', num2str(fovIDs(k))];
    fovMask = (imagePoints(:, 6) == fovIDs(k));
    allFOV.(matlab.lang.makeValidName(fovName)) = imagePoints(fovMask, :);
end
end

function fittedVector = fit_fret_vectors(points)
% FIT_FRET_VECTORS Fits 3D vectors to FRET data using PCA
% Input:  points - Nx5 matrix [x, y, DD_int, DA_int, AA_int]
% Output: fittedVector - 1x3 vector of principal components
%
% Methodology:
% 1. Z-score normalization for outlier detection
% 2. PCA for dimensionality reduction
% 3. Returns first principal component as representative vector[1](@ref)

% Extract intensity data
intensityData = points(:, 3:5);

% Outlier removal using z-scores
zScores = zscore(intensityData);
outlierMask = any(abs(zScores) > 2, 2);  % Threshold = 2 SD
cleanData = intensityData(~outlierMask, :);

% Principal Component Analysis
[coeff, ~, ~] = pca(cleanData);
fittedVector = coeff(:, 1)';  % First principal component
end

function [x_opt, keyParams] = calculate_system_parameters(I_DD, I_DA, I_AA, E_D, E_A)
% CALCULATE_SYSTEM_PARAMETERS Computes FRET calibration parameters
% Inputs: 
%   I_DD - Donor intensities (1x3 vector)
%   I_DA - FRET intensities (1x3 vector)
%   I_AA - Acceptor intensities (1x3 vector)
%   E_D  - Donor FRET efficiencies (1x3 vector)
%   E_A  - Acceptor FRET efficiencies (1x3 vector)
%
% Outputs:
%   x_opt      - Optimization solution vector (9x1)
%   keyParams  - Structure with system parameters (G, a, d, k, eA_eD)
%
% Algorithm:
% 1. Constructs linear system from FRET equations
% 2. Two-stage optimization with constraints
% 3. Physical parameter extraction

% --- Stage 1: Construct Linear System ---
M1 = zeros(3, 3);
M2 = zeros(3, 3);
M3 = zeros(3, 3);

for i = 1:3
    % Equation matrices
    M1(i, :) = [I_DA(i)/I_DD(i), -I_AA(i)/I_DD(i), -1];
    M2(i, :) = [I_DA(i)/I_AA(i), -I_DD(i)/I_AA(i), -1];
    M3(i, :) = [I_DA(i), -I_AA(i), -I_DD(i)];
end

% Compose block-diagonal matrix
A = blkdiag(M1, M1, M2);  % 9x9 system matrix

% Construct right-hand side vector
B = [E_D./(1 - E_D), zeros(1, 3), E_A]';  % 9x1 vector

% --- Stage 2: Two-Step Optimization ---
% Initial parameter estimates
a0 = 0.252;  % Acceptor crosstalk initial
d0 = 0.518;  % Donor crosstalk initial

% First optimization (unconstrained)
x0 = [0.01; 0.01; 0.01];  % Initial guess [G, k, eps_AD]
lb = [0.01; 0.01; 0.01];  % Lower bounds
ub = [ inf;  inf;  inf];    % Upper bounds
options = optimset('Algorithm', 'interior-point', 'Display', 'off');
x1 = fmincon(@(x)initial_objective(x, I_DD, I_DA, I_AA, E_D, E_A, a0, d0),...
             x0, [], [], [], [], lb, ub, [], options);

% Second optimization (constrained)
x0_second = [1/x1(1); a0/x1(1); d0/x1(1); 1; x1(1)*x1(2)+a0; d0-x1(1);...
             x1(3)/a0; (d0*x1(3))/a0; x1(3)];

% Define constraints
lb_second = [0.001; 0.001; 0.001; 1;  1;   -Inf; 0.001; 0.001; 0.001];
ub_second = [0.999; 0.999; 0.999; 1; 10; -0.001;     1;     1;     1];
A_ineq = [ 1, 0, 0, 0, 0, 0, 0, 0, 0;  % x(1) < 1
           0, 1, 0, 0, 0, 0, 0, 0, 0;  % x(2) < 1
           0, 0, 1, 0, 0, 0, 0, 0, 0;  % x(3) < 1
          -1, 1, 0, 0, 0, 0, 0, 0, 0;  % a < 1
          -1, 0, 1, 0, 0, 0, 0, 0, 0;  % d < 1
           0, 0, 0, 0, -1, 0, 0, 0, 0; % x(5) ≥ 1
           0, 0, 0, 0, 0,  1, 0, 0, 0; % x(6) < 0
           0, 0, 0, 0, 0,  0, 1, 0, 0; % x(7) < 1
           0, 0, 0, 0, 0,  0, 0, 1, 0; % x(8) < 1
           0, 0, 0, 0, 0,  0, 0, 0, 1; % x(9) < 1
           ];     
b_ineq = [-1; -1; -1; -eps; -eps;-1; -eps; -1; -1; -1];

A_eq = [0, 0, 0, 1, 0, 0, 0, 0, 0]; 
b_eq = 1;

% Perform constrained optimization
options = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter',...
                       'MaxIterations', 1000, 'StepTolerance', 1e-6,...
                       'ConstraintTolerance', 1e-6);
x_opt = fmincon(@(x)norm(A*x - B), x0_second, A_ineq, b_ineq, A_eq, b_eq,...
                lb_second, ub_second, @nonlincon, options);

% --- Stage 3: Parameter Extraction ---
keyParams.G = 1/x_opt(1);
keyParams.a = keyParams.G * x_opt(2);
keyParams.d = keyParams.G * x_opt(3);
keyParams.k = (x_opt(5) - keyParams.a)*x_opt(1);
keyParams.eA_eD = x_opt(9);
keyParams.Gini = x1(1);
keyParams.kini = x1(2);
keyParams.eps_ADini = x1(3);
end

function f = initial_objective(x, I_DD, I_DA, I_AA, E_D, E_A, a, d)
% INITIAL_OBJECTIVE Objective function for first-stage optimization
% Computes residual sum of squares for initial parameter estimation

G = x(1); k = x(2); eps_AD = x(3);
residuals = zeros(9,1);

    for i = 1:3
        % Residual equations
        residuals(3*i-2) = (G - d)*I_DD(i) + I_DA(i) - (G*k + a)*I_AA(i);
        residuals(3*i-1) = (G*E_D(i)/(E_D(i)-1) - d)*I_DD(i) + I_DA(i) - a*I_AA(i);
        residuals(3*i) = -d*I_DD(i) + I_DA(i) - a*(E_A(i)/eps_AD + 1)*I_AA(i);
    end

f = sum(residuals.^2);  % Sum of squared residuals

end

% Nonlinear constraints
function [c, ceq] = nonlincon(x)            
        c = [];%zeros(4, 1);
        %c(1) = x(3)/x(1) - 1;              % x(8)/x(7) ≤ d' → x(8)/x(7) - d' ≤ 0
        %c(2) = - x(3)/x(1);
        %c(3) = x(2)/x(1) - 1;              % x(2)/x(1) < 1 → x(2)/x(1) - 1 < 0
        %c(4) = 0 - x(2)/x(1);       
        
        ceq = zeros(4, 1);
        ceq(1) = x(8)/x(7) - x(3)/x(1);    % x(8)/x(7) = x(3)/x(1)
        ceq(2) = x(9)/x(7) - x(2)/x(1);    % x(9)/x(7) = x(2)/x(1)
        ceq(3) = -1/x(1)+x(3)/x(1)-x(6);   % 1/x(1) - x(3)/x(1) + x(6) = 0
        ceq(4) = x(8)/x(9) - x(3)/x(2);    % x(8)/x(9) = x(3)/x(2)

end
%%% ==================== END FUNCTIONS ====================== %%%