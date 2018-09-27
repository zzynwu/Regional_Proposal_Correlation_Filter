%   This function runs the BACF tracker on the video specified in "seq".
%   This function borrowed from SRDCF paper. 
%   details of some parameters are not presented in the paper, you can
%   refer to SRDCF/CCOT paper for more details.

function results = run_RPCF(seq, res_path, bSaveImage)

    %   HOG feature parameters
    hog_params.nDim   = 31;
    hog_params.nOrients = 9;
    %   CN feature parameters
    cn_params.tablename = 'CNnorm';
    cn_params.cell_size = 4;
    cn_params.nDim = 10;
    cn_params.useForGray = false;
    %   Grayscale feature parameters
    grayscale_params.colorspace='gray';
    grayscale_params.nDim = 1;
    %   Gray feature parameters
    ic_params.tablename = 'intensityChannelNorm6';
    ic_params.useForColor = false;
    ic_params.cell_size = 4;
    ic_params.nDim = 5;
    %   Global feature parameters 
    params.t_features = {
        ...struct('getFeature',@get_colorspace, 'fparams',grayscale_params),...  % Grayscale is not used as default
        struct('getFeature',@get_fhog,'fparams',hog_params),...
        ...struct('getFeature',@get_table_feature,'fparams',cn_params),...
        ...struct('getFeature',@get_table_feature,'fparams',ic_params),...
    };
    params.t_global.cell_size = 4;                  % Feature cell size
    params.t_global.cell_selection_thresh = 0.75^2; % Threshold for reducing the cell size in low-resolution cases

    %   Search region + extended background parameters
    params.search_area_shape = 'square';    %'proportional';    % the shape of the training/detection window: 'proportional', 'square' or 'fix_padding'
    params.search_area_scale = 5;%5;           % the size of the training/detection area proportional to the target size
    params.filter_max_area   = 50^2;        % the size of the training/detection area in feature grid cells

    %   Color Histogram parameters
    params.n_bins = 32;
    params.learning_rate_pwp = 0.04;
    params.inner_padding = 0.1;  % 0.2
    params.fixed_area = 150^2;

    %   Region proposal parameters
    params.thresh = 0.6;
    params.dist = 25;

    %   Learning parameters
    params.learning_rate       = 0.013;        % learning rate
    params.output_sigma_factor = 1/16;		% standard deviation of the desired correlation output (proportional to target)

    %   Detection parameters
    params.interpolate_response  = 1;        % correlation score interpolation strategy: 0 - off, 1 - feature grid, 2 - pixel grid, 4 - Newton's method
    params.newton_iterations     = 5; %50          % number of Newton's iteration to maximize the detection scores
                    % the weight of the standard (uniform) regularization, only used when params.use_reg_window == 0
    %   Scale parameters
    params.number_of_scales =  3; %5 3 
    params.scale_step       = 1.02; %1.01 1.02

    %   size, position, frames initialization
    params.video_path = seq.path;
    params.wsize    = [seq.init_rect(1,4), seq.init_rect(1,3)];
    params.init_pos = [seq.init_rect(1,2), seq.init_rect(1,1)] + floor(params.wsize/2);
    params.s_frames = seq.s_frames;
    params.no_fram  = seq.len;
    params.seq_st_frame = seq.startFrame;
    params.seq_en_frame = seq.endFrame;
    % params.last_frame = seq.last_frame;

    %   ADMM parameters, # of iteration, and lambda- mu and betha are set in
    %   the main function.
    params.admm_iterations = 2;
    params.admm_lambda = 0.01;

    %   Debug and visualization
    params.visualization = 0;
    params.debug = 0;

    %   Run the main function
    res = RPCF_tracker(params);
    results.type = res.type;
    results.res = res.res;
    results.fps = res.fps;
    disp(['fps: ' num2str(results.fps)])
end