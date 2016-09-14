function [PhaseMasks, TransformedSLMTargets] = SLMPhaseMaskMakerCUDA3D(varargin)
% Lloyd Russell 2015
% Implements a weighted GS algorithm using Martin Persson's HOTlab DLL
% (https://github.com/MartinPersson/HOTlab)

% Parameters
% =========================================================================
% None         : open a load file GUI
% TargetsImage : 512*512 array or { arrays }
% FileName     : string or { strings }
% Points       : n*4 array (x,y,z,I) or { arrays }
% -------------------------------------------------------------------------
% Save         : boolean (default=true)
% SaveName     : string or { strings }, NB must include .tiff extension
% -------------------------------------------------------------------------
% DoTransform  : boolean (default=true)
% TformPath    : string
% -------------------------------------------------------------------------
% Iterations   : integer
% -------------------------------------------------------------------------
% SliceSpacing : integer, step size between slices of input 3D stack
% FocalSlice   : integer, which slice of 3D stack is the normal focal plane
% ManualZ      : integer or 1d array of length n spots or { array }

% Returns
% =========================================================================
% PhaseMasks
% TransformedSLMTargets
% -------------------------------------------------------------------------
% Saved files:
% |- PhaseMasks                     < results folder
% |----- (filename)_CUDAphase.tiff  < the computed phase mask
% |----- InputTargets               < folder of input targets
% |----- TransformedTargets         < folder of transformed targets

% NB: -um puts spots deeper in sample than native focal plan of objective

p = inputParser;
p.addParameter('TargetsImage', []);
p.addParameter('SaveName', '');
p.addParameter('FileName', '');
p.addParameter('Do2DTransform', true);
p.addParameter('Do3DTransform', true);
p.addParameter('TformPath', '20151209_tform_2');
p.addParameter('Save', true);
p.addParameter('Iterations', 100);
p.addParameter('Is3D', false);
p.addParameter('SliceSpacing', 1);
p.addParameter('FocalSlice', 1);
p.addParameter('ManualZ', []);
p.addParameter('Points',[]);
parse(p, varargin{:});

FocalSlice{1} = p.Results.FocalSlice;
SliceSpacing = p.Results.SliceSpacing;

% INPUT: Manual points
% -------------------------------------------------------------------------
if ~any(strcmpi(p.UsingDefaults, 'Points'))
    if any(strcmpi(p.UsingDefaults, 'SaveName'))
        warning('Must provide save name when passing in raw points')
        return
    end
    
    if iscell(p.Results.Points)
        Points = p.Results.Points;
        NumFiles = length(p.Results.Points);
        SaveName = p.Results.SaveName;
    else
        NumFiles = 1;
        Points{1} = p.Results.Points;
        SaveName{1} = p.Results.SaveName;
    end
    
    for idx = 1:NumFiles
        num_targets = numel(Points{idx}(:,1));
        num_planes = range(unique(Points{idx}(:,3)))+1;
        InputTargets{idx} = zeros(512,512,num_planes);
        for i = 1:num_targets
            row = Points{idx}(i,2);
            col = Points{idx}(i,1);
            slice = Points{idx}(i,3) - min(Points{idx}(:,3)) + 1;
            val = Points{idx}(i,4);
            InputTargets{idx}(row,col,slice) = val;
        end
        FocalSlice{idx} = -min(Points{idx}(:,3)) + 1;
    end
    
    Do2DTransform = p.Results.Do2DTransform;
    
    % INPUT: Image array provided
    % -------------------------------------------------------------------------
elseif ~any(strcmpi(p.UsingDefaults, 'TargetsImage'))
    if ~any(strcmpi(p.UsingDefaults, 'SaveName'))  % save name must also be provided
        if iscell(p.Results.TargetsImage)  % passed multiple images
            InputTargets = p.Results.TargetsImage;
            SaveName = p.Results.SaveName;
        else
            NumFiles = 1;
            InputTargets{1} = p.Results.TargetsImage;
            SaveName{1} = p.Results.SaveName;
        end
    else
        warning('Must provide save name when passing in image data')
    end
    Do2DTransform = p.Results.Do2DTransform;
    
    % INPUT: Filepath provided
    % -------------------------------------------------------------------------
elseif ~any(strcmpi(p.UsingDefaults, 'FileName'))
    FileName = p.Results.FileName;
    if iscell(FileName)  % passed multiple filenames
        NumFiles = length(FileName);
        InputTargets = cell(NumFiles,1);
        
        for idx = 1:length(FileName)
            info = imfinfo(FileName{idx});
            num_images = numel(info);
            InputTargets{idx} = zeros(512,512,num_images);
            for k = 1:num_images
                InputTargets{idx}(:,:,k) = imread(FileName{idx}, k);
            end
        end
    else
        NumFiles = 1;
        info = imfinfo(FileName);
        num_images = numel(info);
        InputTargets{1} = zeros(512,512,num_images);
        for k = 1:num_images
            InputTargets{1}(:,:,k) = imread(FileName, k);
        end
        FileName = {FileName};
    end
    SaveName = FileName;
    Do2DTransform = p.Results.Do2DTransform;
    
    
    % INPUT: none, open load file dialog
    % -------------------------------------------------------------------------
elseif any(strcmpi(p.UsingDefaults, 'TargetsImage')) && any(strcmpi(p.UsingDefaults, 'FileName'))
    [FileName, PathName] = uigetfile('*.tif*', 'Select the targets file(s)', 'MultiSelect', 'on');
    if PathName == 0
        warning('Cancelled')
        return
    end
    transform_choice = questdlg('Transform points?', 'User input required', 'Yes', 'No', 'Yes');
    if strcmpi(transform_choice, 'Yes')
        Do2DTransform = true;
    else
        Do2DTransform = false;
    end
    
    % load images into array
    cd(PathName);
    InputTargets = {};
    if iscell(FileName)  % selected multiple files
        NumFiles = length(FileName);
        InputTargets = cell(NumFiles,1);
        for idx = 1:length(FileName)
            info = imfinfo(FileName{idx});
            num_images = numel(info);
            InputTargets{idx} = zeros(512,512,num_images);
            for k = 1:num_images
                InputTargets{idx}(:,:,k) = imread(FileName{idx}, k);
            end
        end
    else
        NumFiles = 1;
        info = imfinfo(FileName);
        num_images = numel(info);
        InputTargets{1} = zeros(512,512,num_images);
        for k = 1:num_images
            InputTargets{1}(:,:,k) = imread(FileName, k);
        end
        FileName = {FileName};
    end
    SaveName = FileName;
end


Do3DTransform = p.Results.Do3DTransform;


% rebuild input targets if provided manual Z coordinates
if ~any(strcmpi(p.UsingDefaults, 'ManualZ'))
    if iscell(p.Results.ManualZ);
        manualZ = p.Results.ManualZ;
    else
        manualZ{1} =  p.Results.ManualZ;
    end
    if numel(manualZ) < NumFiles
        manualZ = repmat(manualZ,1,NumFiles);
    end
    
    for idx = 1:NumFiles
        % Find target coordinates
        [y, x, z] = ind2sub(size(InputTargets{idx}),find(InputTargets{idx}));
        
        % get target intensities
        [~,~,I] = find(InputTargets{idx});
        
        % replace Z with manual Z
        if numel(manualZ{idx}) ~= numel(z)
            manualZ{idx} = repmat(manualZ{idx},1,numel(z));
        end
        z = manualZ{idx};
        
        % build new input image
        num_targets = numel(x);
        num_planes = range(unique(z))+1;
        InputTargets{idx} = zeros(512,512,num_planes);
        for i = 1:num_targets
            row = y(i);
            col = x(i);
            slice = z(i) - min(z) + 1;
            val = I(i);
            InputTargets{idx}(row,col,slice) = val;
        end
        FocalSlice{idx} = -min(z) + 1;
    end
end

% Make output directories
if p.Results.Save
    if ~exist('PhaseMasks', 'dir')
        mkdir('PhaseMasks');
    end
    if ~exist(['PhaseMasks' filesep 'InputTargets'], 'dir')
        mkdir(['PhaseMasks' filesep 'InputTargets']);
    end
    if ~exist(['PhaseMasks' filesep 'TransformedTargets'], 'dir')
        mkdir(['PhaseMasks' filesep 'TransformedTargets']);
    end
end

% Load previously calibrated transform
if Do2DTransform
    load(p.Results.TformPath)
end

% Load HOTlab DLL
if ~libisloaded('GenerateHologramCUDA')
    loadlibrary('GenerateHologramCUDA')
    %     libfunctions('GenerateHologramCUDA')
end

% 'StartCUDAandSLM' parameters
EnableSLM   = 0;
TrueFrames  = 6;
deviceId    = 0;
h_pSLMstart = (rand(512, 512)-0.5)*2*pi;  % the starting phase mask
LUT         = [];

% Start CUDA (and SLM)
[cuda_error1,~,~] = calllib('GenerateHologramCUDA','startCUDAandSLM',...
    EnableSLM, h_pSLMstart, LUT, TrueFrames, deviceId);

% Process all files (do transforms, make phase masks and save images)
PhaseMasks = cell(NumFiles,1);
TransformedSLMTargets = cell(NumFiles,1);
for f = 1:NumFiles
    tic;
    % Find target coordinates
    [y_i, x_i, z_i] = ind2sub(size(InputTargets{f}),find(InputTargets{f}));
    
    % get target intensities
    [~,~,I] = find(InputTargets{f});
    
    % adjust z coordinates
    z_i = (z_i - FocalSlice{f}) * SliceSpacing;
    
    % get how many different planes
    im_dims = size(InputTargets{f});
    num_dims = length(im_dims);
    if num_dims > 2
        num_planes_input = im_dims(3);
    else
        num_planes_input = 1;
    end
    
    % % for testing
    % figure; hold on; axis square;
    % plot3(x_i, y_i, z_i, 'ro')
    % xlim([0 511]); ylim(xlim)
    
    if Do2DTransform
        % Transform target XY coordinates from 2P to SLM space
        [x, y] = transformPointsForward(tform, x_i, y_i);
    else
        x = x_i; y = y_i;
    end
    
    if Do3DTransform
        % NOTE, TODO: all the below should be replaced with a simple application
        % of a precomputed transformation. See for ideas:
        %    http://uk.mathworks.com/help/vision/examples/3-d-point-cloud-registration-and-stitching.html
        % But note rigid body only at the moment (not suitable for this)
        % Just compute the transform matrix manually.
        
        
        
        
        x = x_i;
        y = y_i;

        % FOR XY ADJUSTMENTS ACTUAL Z UM ARE NOT YET REQUIRED SO LONG AS
        % CALIBRATIONS ARE DONE IN A.U
        
        
        % DO XY OFFSET WITH Z 
        % ------======-------
        % This is performed in 2P imaging space. SO the offset in X and Y
        % must be measured in these units
        
        % calculate offset in imaging pixels - doesn't need to be this
        % complicated - but I imaged at 8x to get higher resolution...
        z_au = [100 0 -100];
        x_offset_um = [108*0.13 0 -102*0.13];
        y_offset_um = [-11*0.13 0 11*0.13];
        working_fov = 399.5;  % and 512 pixels
        um_per_pixel = working_fov / 512;
        x_offset_pixels = x_offset_um * um_per_pixel;
        y_offset_pixels = y_offset_um * um_per_pixel;
       
        x_offsets = interp1(z_au, x_offset_pixels, z_i, 'spline');
        y_offsets = interp1(z_au, y_offset_pixels, z_i, 'spline');
                
        x = x - x_offsets;
        y = y - y_offsets;
        
        
        % DO XY SCALING WITH Z  (should this be before or after the offset?)
        % ------=======-------
        sample_points = [100 50 0 -50 -100];  % au
        measured_scales = [162.56 159.3 155.4 150.3 143.3];
        zero_index = find(sample_points==0);
        norm_sizes = measured_scales / measured_scales(zero_index);
        z_scale_factors = interp1(sample_points, norm_sizes, z_i, 'spline');
        
        x = ((x-256) ./ z_scale_factors) + 256;
        y = ((y-256) ./ z_scale_factors) + 256;


        % AFTER PERFORMING THE XY SCALE AND OFFSET CORRECTIONS IN 2P
        % IMAGING SPACE, TRANSFORM INTO SLM PIXELS
        %----------------=========----------------------------------
        if Do2DTransform
            [x, y] = transformPointsForward(tform, x, y);
        end
        
        % XY COORDINATES ARE NOW FINALISED SO CAN DO Z...
        
        
        % DO Z AU TO UM CONVERSION
        % -----========-----------
        % calculate radial distance from zero order
        dist = sqrt(((x-256).^2) + ((y-256).^2));
        
        AU = [100 50 0 -50 -100];  % au
        DIST = [17.248 57.5305 97.8407];  % radial dist (SLM space)
        UM = [120 122 128
              50 58 62
              8 4 2
              -64 -62 -62
              -116 -118 -116];  % measured
        
        AU_RANGE = AU;
        INPUT_DIST = dist;
        EST_UMs = interp2(DIST,AU,UM, INPUT_DIST, AU_RANGE, 'spline');
        
        % do each spot in turn becasue xy has an efefct
        required_au = [];
        for i = 1:numel(x)
            desired_offset = z_i(i);
            required_au(i) = interp1(EST_UMs(:,i),AU_RANGE,desired_offset, 'spline');
        end
        z = required_au;
        z(z_i==0 & dist<17) = 0;  % prevent offset for intended 0 and close spots
        
    else
        z = z_i;
    end
    
    
    % round
    x = round(x);
    y = round(y);
    z = round(z);  % neccessary to round z units?
    
    
    % Check for and remove targets outside of SLM space
    if any(x<1 | x>512 | y<1 | y>512)
        warning('At least one transformed target spot is outside of addressable SLM space. The offending spots will be removed.')
        idx_remove = (x<1 | x>512 | y<1 | y>512);
        x(idx_remove) = [];
        y(idx_remove) = [];
        z(idx_remove) = [];
        I(idx_remove) = [];
    end
    
    % get how many different planes
    num_planes_output = range(unique(z))+1;
    
    % Build transformed targets image
    TransformedSLMTargets{f} = uint8(zeros(512, 512, num_planes_output));
    for i = 1:length(x)
        TransformedSLMTargets{f}(y(i), x(i), z(i)-min(z)+1) = I(i);
    end
    
    % 'GenerateHologramCUDA' parameters
    h_test         = [];
    h_pSLM         = zeros(512, 512);
    x_spots        = x - 256;  % subtract 256 because centre of image is 0,0 (not 256,256)
    y_spots        = y - 256;  % subtract 256 because centre of image is 0,0 (not 256,256)
    z_spots        = z;
    I_spots        = I;  % intensity comes from input image
    N_spots        = length(x_spots);
    N_iterations   = p.Results.Iterations;
    h_obtainedAmps = [];
    if any(z_spots ~= 0)
        method = 1;
    else
        method = 2;
    end
    % 0: Complex addition of "Lenses and Prisms", no optimization (3D)
    % 1: Weighted Gerchberg-Saxton algorithm using Fresnel propagation (3D)
    % 2: Weighted Gerchberg-Saxton algorithm using fast fourier transforms (2D)
    
    % Apply corrections
    % [cuda_error, h_AberrationCorr, h_LUTPolCoeff, h_LUT_uc] = calllib('GenerateHologramCUDA', 'Corrections',...
    %     UseAberrationCorr, h_AberrationCorr, UseLUTPol, PolOrder, h_LUTPolCoeff, saveAmplitudes, alpha, DCborderWidth, UseLUT, h_LUT_uc);
    
    % Generate Hologram
    [cuda_error2,~,h_pSLM,~,~,~,~,h_obtainedAmps] = calllib('GenerateHologramCUDA','GenerateHologram',...
        h_test, h_pSLM, x_spots, y_spots, z_spots, I_spots, N_spots, N_iterations, h_obtainedAmps, method);
    
    % Convert from 8bit to 16bit
    PhaseMasks{f} = uint16(double(h_pSLM)/255*65535);
    
    % Save
    if p.Results.Save
        imwrite(PhaseMasks{f}, ['PhaseMasks' filesep strrep(SaveName{f}, '.tif', '_CUDAphase.tiff')]);
        for k = 1:num_planes_input
            if k == 1;
                writemode = 'overwrite';
            else
                writemode = 'append';
            end
            imwrite(InputTargets{f}(:,:,k), ['PhaseMasks' filesep 'InputTargets' filesep strrep(SaveName{f}, '.tif', '_InputTargets.tif')], 'writemode', writemode);
        end
        
        for k = 1:num_planes_output
            if k == 1;
                writemode = 'overwrite';
            else
                writemode = 'append';
            end
            imwrite(TransformedSLMTargets{f}(:,:,k), ['PhaseMasks' filesep 'TransformedTargets' filesep strrep(SaveName{f}, '.tif', '_TransformedTargets.tif')], 'writemode', writemode);
        end
    end
    time_stop = toc;
    disp(['[' num2str(f) ' of ' num2str(NumFiles) '] Done in ' num2str(round(time_stop*1000*1000)/1000) ' ms'])
end

% Stop CUDA (and SLM)
[cuda_error3] = calllib('GenerateHologramCUDA','stopCUDAandSLM');
% unloadlibrary('GenerateHologramCUDA');  % tidy up, unload dll
