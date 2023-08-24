%% Modified by Mahe Chen, Alessandro Luchetti and Andrew Mocle
%% Clear Workspace
clear; clc; close all;
%% Select Paths
addpath(fullfile(pwd, 'MatlabPath_CNMFe_and_NoRMCorre'));
CellReg_directory=fileparts(which('CellReg_setup.m'));
addpath(genpath(fullfile(CellReg_directory,'CellReg')));
addpath(genpath(fullfile(CellReg_directory,'GUI')));
addpath(genpath(fullfile(CellReg_directory,'Helper')));
savepath();
%% Settings
settings.spatial_downsampling = 2; % Recommended range: 2 - 4. Downsampling significantly increases computational speed, but verify it does not
settings.path = uigetdir;
settings.equipment = 'V4'; % Equipment used for imaging.
settings.Sessions = dir(settings.path);
settings.Sessions = settings.Sessions(3:end, :);
settings.analysis_time = 'Result';
isnonrigid = true; % If true, performs non-rigid registration (slower). If false, rigid alignment (faster). Non-rigid is preferred within sessions.
plotFlag = false; % Plot the results of motion correction
ROIflag = false; % Choose true if you want to select a specific ROI in the FOV for each separate session. Pixels outside of the FOV will be deemed zero.
replaceRGBVideo = false; % Choose true if you want to replace RGB videos by their grayscale version
checkMotCorr = false;
nSessions = size(settings.Sessions, 1);
save(fullfile(settings.path, 'settings.mat'), 'settings', '-v7.3')
%% Video Processing
disp('Applying motion correction on single sessions.');
for i = 1:nSessions
    %% Apply motion correction and concatenate videos
    sessionPath = fullfile(settings.path, settings.Sessions(i).name);
    cd(sessionPath)
    ms = msGenerateVideoObjConcat(sessionPath, settings.equipment, replaceRGBVideo, 'msCam');
    ms.FrameRate = round(1/(median(diff(ms.time), "omitmissing") / 1000));
    ms.equipment = settings.equipment;
    if settings.FrameRate ~= ms.FrameRate
        disp(['Difference in Frame Rates Detected!']);
    end
    settings.FrameRate = ms.FrameRate;
    settings.normalize = true;
    settings.align = false;
    ms.analysis_time = settings.analysis_time;
    ms.ds = settings.spatial_downsampling;
    analysisFolder = fullfile(sessionPath, settings.analysis_time);
    mkdir(analysisFolder);
    save(fullfile(ms.dirName, 'ms.mat'), 'ms');
    disp(['Working on Session: ' num2str(i) ' of ' num2str(nSessions)])
    ms = msNormCorreConcat(ms, isnonrigid, ROIflag, plotFlag, checkMotCorr);
    save(fullfile(ms.dirName, 'ms.mat'), 'ms');
    clear ms
    %% Video Normalization
    if settings.normalize
        inputVideo = VideoReader(fullfile(settings.path, settings.Sessions(i).name, settings.analysis_time, 'msvideo.avi'));
        baselineFrame = zeros(inputVideo.Height, inputVideo.Width);
        frameCount = 0;
        while hasFrame(inputVideo)
            frame = readFrame(inputVideo);
            frame = im2double(frame);
            baselineFrame = baselineFrame + frame;
            frameCount = frameCount + 1;
        end
        baselineFrame = baselineFrame / frameCount;
        inputVideo = VideoReader(fullfile(settings.path, settings.Sessions(i).name, settings.analysis_time, 'msvideo.avi'));
        outputVideo = VideoWriter(fullfile(settings.path, settings.Sessions(i).name, settings.analysis_time, 'normvideo.avi'));
        outputVideo.FrameRate = inputVideo.FrameRate;
        open(outputVideo);
        while hasFrame(inputVideo)
            frame = readFrame(inputVideo);
            frame = im2double(frame);
            frame = imfilter(frame, fspecial('unsharp'));
            frame = imgaussfilt(frame, 2);
            deltaF = frame - baselineFrame;
            deltaF = denoiseTV(deltaF, 0.1, 10);
            deltaF = adapthisteq(deltaF);
            deltaF = max(0, min(deltaF, 1)); 
            deltaF = uint8(deltaF * 255);
            writeVideo(outputVideo, deltaF);
        end
        close(outputVideo);
    end
end
%% CNMF-E
for j = 1:1
    %% choose data
    neuron = Sources2D();
    neuron.select_data(fullfile(settings.path, settings.Sessions(j).name, settings.analysis_time, 'msvideo.avi'));  %if nam is [], then select data interactively
    savename = "data.mat";

    %% parameters
    %min_pnr and min_corr (res as well) typicall 5, 0.5
    % -------------------------    COMPUTATION    -------------------------  %
    pars_envs = struct('memory_size_to_use', 4, ...   % GB, memory space you allow to use in MATLAB
        'memory_size_per_patch', 2, ...   % GB, space for loading data within one patch
        'patch_dims', [64, 64]);  %GB, patch size

    % -------------------------      SPATIAL      -------------------------  %
    gSig = 3;           % pixel, gaussian width of a gaussian kernel for filtering the data. 0 means no filtering
    gSiz = 15;          % pixel, neuron diameter
    ssub = settings.spatial_downsampling;           % spatial downsampling factor
    with_dendrites = false;   % with dendrites or not
    if with_dendrites
        % determine the search locations by dilating the current neuron shapes
        updateA_search_method = 'dilate';  %#ok<UNRCH>
        updateA_bSiz = 5;
        updateA_dist = neuron.options.dist;
    else
        % determine the search locations by selecting a round area
        updateA_search_method = 'ellipse';
        updateA_dist = 5;
        updateA_bSiz = neuron.options.dist;
    end
    spatial_constraints = struct('connected', true, 'circular', false);  % you can include following constraints: 'circular'
    spatial_algorithm = 'hals';
    
    % -------------------------      TEMPORAL     -------------------------  %
    load(fullfile(settings.path, settings.Sessions(j).name, 'ms.mat'))
    Fs = ms.FrameRate;             % frame rate
    clear ms
    tsub = 5;           % temporal downsampling factor
    deconv_flag = true;     % run deconvolution or not 
    deconv_options = struct('type', 'ar1', ... % model of the calcium traces. {'ar1', 'ar2'}
        'method', 'foopsi', ... % method for running deconvolution {'foopsi', 'constrained', 'thresholded'}
        'smin', -3, ...         % minimum spike size. When the value is negative, the actual threshold is abs(smin)*noise level
        'optimize_pars', true, ...  % optimize AR coefficients
        'optimize_b', true, ...% optimize the baseline);
        'max_tau', 100);    % maximum decay time (unit: frame);

    nk = 3;             % detrending the slow fluctuation. usually 1 is fine (no detrending)
    % when changed, try some integers smaller than total_frame/(Fs*30)
    detrend_method = 'spline';  % compute the local minimum as an estimation of trend.

    % -------------------------     BACKGROUND    -------------------------  %
    bg_model = 'ring';  % model of the background {'ring', 'svd'(default), 'nmf'}
    nb = 1;             % number of background sources for each patch (only be used in SVD and NMF model)
    ring_radius = 20;  % when the ring model used, it is the radius of the ring used in the background model.
    %otherwise, it's just the width of the overlapping area
    num_neighbors = []; % number of neighbors for each neuron
    %bg_ssub = 1;        % downsample background for a faster speed 

    % -------------------------      MERGING      -------------------------  %
    show_merge = false;  % if true, manually verify the merging step
    merge_thr = 0.65;     % thresholds for merging neurons; [spatial overlap ratio, temporal correlation of calcium traces, spike correlation]
    method_dist = 'max';   % method for computing neuron distances {'mean', 'max'}
    dmin = 5;       % minimum distances between two neurons. it is used together with merge_thr
    dmin_only = 2;  % merge neurons if their distances are smaller than dmin_only.
    merge_thr_spatial = [0.8, 0.4, -inf];  % merge components with highly correlated spatial shapes (corr=0.8) and small temporal correlations (corr=0.1)

    % -------------------------  INITIALIZATION   -------------------------  %
    K = [];             % maximum number of neurons per patch. when K=[], take as many as possible.
    min_corr = 0.4;     % minimum local correlation for a seeding pixel
    min_pnr = 4;       % minimum peak-to-noise ratio for a seeding pixel
    min_pixel = gSig^2;      % minimum number of nonzero pixels for each neuron
    bd = 5;             % number of rows/columns to be ignored in the boundary (mainly for motion corrected data)l
    frame_range = [];   % when [], uses all frames
    save_initialization = false;    % save the initialization procedure as a video.
    use_parallel = true;    % use parallel computation for parallel computing
    show_init = false;   % show initialization results
    choose_params = false; % manually choose parameters
    center_psf = true;  % set the value as true when the background fluctuation is large (usually 1p data)
    % set the value as false when the background fluctuation is small (2p)

    % -------------------------  Residual   -------------------------  %
    min_corr_res = 0.6; %6
    min_pnr_res = 7; %7
    seed_method_res = 'auto';  % method for initializing neurons from the residual
    update_sn = true;

    % ----------------------  WITH MANUAL INTERVENTION  --------------------  %
    with_manual_intervention = true;

    % -------------------------  FINAL RESULTS   -------------------------  %
    save_demixed = true;    % save the demixed file or not
    kt = 3;                 % frame intervals

    % -------------------------    UPDATE ALL    -------------------------  %
    neuron.updateParams('gSig', gSig, ...       % -------- spatial --------
        'gSiz', gSiz, ...
        'ring_radius', ring_radius, ...
        'ssub', ssub, ...
        'search_method', updateA_search_method, ...
        'bSiz', updateA_bSiz, ...
        'dist', updateA_bSiz, ...
        'spatial_constraints', spatial_constraints, ...
        'spatial_algorithm', spatial_algorithm, ...
        'tsub', tsub, ...                       % -------- temporal --------
        'deconv_flag', deconv_flag, ...
        'deconv_options', deconv_options, ...
        'nk', nk, ...localBG
        'detrend_method', detrend_method, ...
        'background_model', bg_model, ...       % -------- background --------
        'nb', nb, ...
        'ring_radius', ring_radius, ...
        'num_neighbors', num_neighbors, ...
        'merge_thr', merge_thr, ...             % -------- merging ---------
        'dmin', dmin, ...
        'method_dist', method_dist, ...
        'min_corr', min_corr, ...               % ----- initialization -----
        'min_pnr', min_pnr, ...
        'min_pixel', min_pixel, ...
        'bd', bd, ...
        'center_psf', center_psf);
    neuron.Fs = Fs;

    %% distribute data and be ready to run source extraction
    neuron.getReady(pars_envs);

    %% initialize neurons from the video data within a selected temporal range
    if choose_params
        % change parameters for optimized initialization
        [gSig, gSiz, ring_radius, min_corr, min_pnr] = neuron.set_parameters();
    end

    [center, Cn, PNR] = neuron.initComponents_parallel(K, frame_range, save_initialization, use_parallel);
    neuron.compactSpatial();
    if show_init
        figure();
        ax_init= axes();
        imagesc(Cn, [0, 1]); colormap gray;
        hold on;
        plot(center(:, 2), center(:, 1), '.r', 'markersize', 10);
    end

    %% estimate the background components
    neuron.update_background_parallel(use_parallel);
    %neuron_init = neuron.copy();

    %% merge neurons and update spatial/temporal components
    neuron.merge_neurons_dist_corr(show_merge);
    neuron.merge_high_corr(show_merge, merge_thr_spatial);

    %% update spatial components

    %% pick neurons from the residual
    [center_res, Cn_res, PNR_res] =neuron.initComponents_residual_parallel([], save_initialization, use_parallel, min_corr_res, min_pnr_res, seed_method_res);
    if show_init
        axes(ax_init);
        plot(center_res(:, 2), center_res(:, 1), '.g', 'markersize', 10);
    end

    %neuron_init_res = neuron.copy();

    %% Initialize AlessandROI and Update Spatial & Temporal Components
    % update spatial
    if update_sn
        neuron.update_spatial_parallel(use_parallel, true);
        update_sn = false;
    else
        neuron.update_spatial_parallel(use_parallel);
    end
    % merge neurons based on correlations 
    neuron.merge_high_corr(show_merge, merge_thr_spatial);
    neuron.update_background_parallel(use_parallel);
    for m=1:3
        neuron.update_spatial_parallel(use_parallel);
        neuron.trimSpatial(0.01, 3); % for each neuron, apply imopen first and then remove pixels that are not connected with the center
        neuron.compactSpatial();    % run this line if neuron shapes are circular 
        % update temporal
        neuron.update_temporal_parallel(use_parallel);
        % delete bad neurons
        neuron.remove_false_positives();
        % merge neurons based on temporal correlation + distances 
        neuron.merge_neurons_dist_corr(show_merge);
        [center_res, Cn_res, PNR_res] =neuron.initComponents_residual_parallel([], save_initialization, use_parallel, min_corr_res, min_pnr_res, seed_method_res);
    end
    %% Run AlessandROI to Remove False Positives
    neuron.options.spatial_algorithm = 'nnls';
    vid_path = fullfile(settings.path, settings.Sessions(j).name, settings.analysis_time);
    if with_manual_intervention
        show_merge = true;
        neuron.orderROIs('snr');   % order neurons in different ways {'snr', 'decay_time', 'mean', 'circularity'}
        neuron.Coor = neuron.get_contours(0.6);
        neuron.AlessandROI(vid_path, [], neuron.C_raw);

        % merge closeby neurons
        neuron.merge_close_neighbors(true, dmin_only);

        % delete neurons
        tags = neuron.tag_neurons_parallel();  % find neurons with fewer nonzero pixels than min_pixel and silent calcium transients
        ids = find(tags>0); 
        if ~isempty(ids)
            neuron.AlessandROI(vid_path, ids, neuron.C_raw);
        end
    end
    neuron.update_spatial_parallel(use_parallel);
    neuron.trimSpatial(0.01, 3); % for each neuron, apply imopen first and then remove pixels that are not connected with the center
    neuron.compactSpatial()
    neuron.update_temporal_parallel(use_parallel);

    K = size(neuron.A,2);
    tags = neuron.tag_neurons_parallel();  % find neurons with fewer nonzero pixels than min_pixel and silent calcium transients
    neuron.remove_false_positives();
    neuron.merge_neurons_dist_corr(show_merge);
    neuron.merge_high_corr(show_merge, merge_thr_spatial);

    if K~=size(neuron.A,2)
        neuron.update_spatial_parallel(use_parallel);
        neuron.trimSpatial(0.01, 3); % for each neuron, apply imopen first and then remove pixels that are not connected with the center
        neuron.compactSpatial()
        neuron.update_temporal_parallel(use_parallel);
        neuron.remove_false_positives();
    end
    %% run more iterations
    %neuron.update_background_parallel(use_parallel);
    neuron.update_spatial_parallel(use_parallel);
    neuron.trimSpatial(0.01, 3); % for each neuron, apply imopen first and then remove pixels that are not connected with the center
    neuron.compactSpatial()
    neuron.update_temporal_parallel(use_parallel);

    K = size(neuron.A,2);
    tags = neuron.tag_neurons_parallel();  % find neurons with fewer nonzero pixels than min_pixel and silent calcium transients
    neuron.remove_false_positives();
    neuron.merge_neurons_dist_corr(show_merge);
    neuron.merge_high_corr(show_merge, merge_thr_spatial);

    if K~=size(neuron.A,2)
        neuron.update_spatial_parallel(use_parallel);
        neuron.trimSpatial(0.01, 3); % for each neuron, apply imopen first and then remove pixels that are not connected with the center
        neuron.compactSpatial()
        neuron.update_temporal_parallel(use_parallel);
        neuron.remove_false_positives();
    end
    
    if with_manual_intervention
        show_merge = true;
        neuron.orderROIs('snr');   % order neurons in different ways {'snr', 'decay_time', 'mean', 'circularity'}
        neuron.AlessandROI(vid_path, [], neuron.C_raw);

        % merge closeby neurons
        neuron.merge_close_neighbors(true, dmin_only);

        % delete neurons
        tags = neuron.tag_neurons_parallel();  % find neurons with fewer nonzero pixels than min_pixel and silent calcium transients
        ids = find(tags>0); 
        if ~isempty(ids)
            neuron.AlessandROI(vid_path, ids, neuron.C_raw);
        end
    end
    %% save the workspace for future analysis
    neuron.orderROIs('snr');
    [spath, ~, ~] = fileparts(neuron.file);
    save_path = strcat(spath, "/", savename);
    st = struct(neuron); save(save_path, 'st');
     %% show neuron contours
    Coor = neuron.show_contours(0.6);
    saveas(gcf, fullfile(settings.path, int2str(j)), 'jpeg');
    %% Prepare Cellreg data
    s = size(neuron.Cn);
    h = s(1); w = s(2);
    s1 = size(neuron.A); n = s1(2);
    d = reshape(full(neuron.A), h, w, n);
    cellreg_data = permute(d, [3,1,2]);
    save(fullfile(settings.path, settings.Sessions(j).name, settings.analysis_time, sprintf('cellreg_data%1i.mat', j)), 'cellreg_data');
end
%% Optional Alignment Step
if settings.align
    p = '/media/jfadmin/LaCie/Mahe';
    asim = '/media/jfadmin/LaCie/Mahe/data';
    files = ["training" "recent" "remote"];
    mice = ["038-01" "038-03" "038-04" "039-02" "039-04"];
    for i=1:5
        animal={};
        for k=1:length(files)
            load(fullfile(p, files(k), int2str(i), 'ms.mat'));
            animal{k} = ms;
        end
        disp(['Aligning training(1) recent(2) remote(3) of' mice(i)]);
        [settings.AllAlignment,settings.AllCorrelation]=AlignAcrossSessions(animal);
        settings.refAverCorr = nanmax(nanmean(settings.AllCorrelation));
        settings.refSession = 1;
        num_videos = 3;
        settings.FinalAlignment = settings.AllAlignment(settings.refSession,:);

        for a=1:3
            for b=1:5
                load(fullfile(asim, mice(b), strcat(files(a), '.mat')))
                dim = size(cellreg_data);
                num_neurons = dim(1);
                for c=1:num_neurons
                    frame = squeeze(cellreg_data(c,:,:));
                    adjusted_frame = imwarp(frame,settings.FinalAlignment{a},'OutputView',imref2d(size(frame)));
                    adjusted_cellreg_data(c,:,:) = adjusted_frame;
                end
                save(fullfile(asim, mice(b), strcat('adjusted_', files(a), '.mat')), 'adjusted_cellreg_data');
                clear adjusted_cellreg_data;
            end 
        end
    end
end