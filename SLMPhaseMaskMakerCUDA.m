function SLMPhaseMaskMakerCUDA(varargin)
% Lloyd Russell 20151028
% Implements a weighted GS algorithm using Martin Perssons HOTlab DLL
% (https://github.com/MartinPersson/HOTlab)
% 
% Input arguments can be:
%   none, an open file dialog box will open
%   a path to an image file on disk
%   loaded image data, save name must be provided as second argument

% --- Load targets image ---
if nargin == 0  % open file dialog if no filepath provided
    [file_name,path_name] = uigetfile('*.tif*', 'Select the targets file');
    file_path = [path_name filesep file_name];
    cd(path_name);
    targets_img = imread(file_path);
    transform_choice = questdlg('Transform points?', 'User input required', 'Yes', 'No', 'Yes');
elseif ischar(varargin{1})  % filepath provided
    file_path = varargin{1};
    targets_img = imread(file_path);
    transform_choice = 'Yes';
else  % assume data array provided first, and savename second
    targets_img = varargin{1};
    file_path = varargin{2};
    transform_choice = 'Yes';
end
if ~exist('PhaseMasks', 'dir')
    mkdir('PhaseMasks');
end
[~, img_name] = fileparts(file_path);

% --- Find target coordinates ---
[y_i, x_i, I] = find(targets_img);

% --- Load previously calibrated transform ---
load('20151008_tform_2P_SLM')

% --- Convert target *coordinates* directly from 2P to SLM space ---
switch transform_choice
    case 'Yes'
        [x, y] = transformPointsForward(tform, x_i, y_i);
        x = round(x); y = round(y);
    case 'No'
        x = x_i; y = y_i;
end

% --- Check for and remove targets outside of SLM space
if any(x<1) || any(x>512) || any(y<1) || any(y>512)
    warning('At least one transformed target spot is outside of addressable SLM space. The offending spots will be removed.')
    idx_remove = (x<1 | x>512 | y<1 | y>512);
    x(idx_remove) = []; y(idx_remove) = [];
end

% --- Build transformed targets image ---
SLM_targets = uint8(zeros(512,512));
for i = 1:length(x)
    SLM_targets(y(i),x(i)) = I(i);
end
imwrite(SLM_targets, [img_name '_Transformed.tif']);

% --- load HOTlab DLL ---
if ~libisloaded('GenerateHologramCUDA')
    loadlibrary('GenerateHologramCUDA')
%     libfunctions('GenerateHologramCUDA')
%     libfunctionsview('GenerateHologramCUDA')
end

% --- start SLM params ---
EnableSLM = 0;
TrueFrames = 6;
deviceId = 0;
h_pSLMstart = (rand(512)-0.5)*2*pi;  % the starting phase mask (seed)
LUT = [];

% --- generate hologram variables ---
h_test = [];
h_pSLM = zeros(512,512);
x_spots = x - 256 - 1;  % subtract 256 because centre of image is 0,0 (not 256,256), subtract 1 because zero indexing
y_spots = y - 256 - 1;  % subtract 256 because centre of image is 0,0 (not 256,256)
z_spots = zeros(length(x_spots),1);  % z planes all set to 0
I_spots = I;  % intensity comes from input image
N_spots = length(x_spots);
N_iterations = 1000;
h_obtainedAmps = [];
method = 2;  % 0: Complex addition of "Lenses and Prisms", no optimization (3D)
             % 1: Weighted Gerchberg-Saxton algorithm using Fresnel propagation (3D)
             % 2: Weighted Gerchberg-Saxton algorithm using fast fourier transforms (2D)

% --- Start CUDA (and SLM) ---
[cuda_error1,~,~] = calllib('GenerateHologramCUDA','startCUDAandSLM',...
    EnableSLM, h_pSLMstart, LUT, TrueFrames, deviceId);

% % --- Apply corrections ---
% [cuda_error, h_AberrationCorr, h_LUTPolCoeff, h_LUT_uc] = calllib('GenerateHologramCUDA', 'Corrections',...
%     UseAberrationCorr, h_AberrationCorr, UseLUTPol, PolOrder, h_LUTPolCoeff, saveAmplitudes, alpha, DCborderWidth, UseLUT, h_LUT_uc);

% --- Generate Hologram ---
[cuda_error2,~,h_pSLM,~,~,~,~,h_obtainedAmps] = calllib('GenerateHologramCUDA','GenerateHologram',...
    h_test, h_pSLM, x_spots, y_spots, z_spots, I_spots, N_spots, N_iterations, h_obtainedAmps, method);

% --- Stop CUDA (and SLM) ---
[cuda_error3] = calllib('GenerateHologramCUDA','stopCUDAandSLM');
unloadlibrary('GenerateHologramCUDA');  % tidy up, unload dll

% --- Convert ---
h_pSLM_16 = uint16(double(h_pSLM)/255*65535);

% --- Save ---
imwrite(h_pSLM_16, ['PhaseMasks' filesep img_name '_CUDAphase.tiff']);
