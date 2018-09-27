% This function implements the BACF tracker.

function [results] = RPCF_tracker(params)

%   Setting parameters for local use.
search_area_scale   = params.search_area_scale;
output_sigma_factor = params.output_sigma_factor;
learning_rate       = params.learning_rate;
filter_max_area     = params.filter_max_area;
nScales             = params.number_of_scales;
scale_step          = params.scale_step;
interpolate_response = params.interpolate_response;

features    = params.t_features;
video_path  = params.video_path;
s_frames    = params.s_frames;
pos         = floor(params.init_pos);
target_sz   = floor(params.wsize);

visualization  = params.visualization;
debug          = params.debug;
num_frames     = params.no_fram;
init_target_sz = target_sz;

%set the feature ratio to the feature-cell size
featureRatio = params.t_global.cell_size;
search_area = prod(init_target_sz / featureRatio * search_area_scale);

% when the number of cells are small, choose a smaller cell size
if isfield(params.t_global, 'cell_selection_thresh')
    if search_area < params.t_global.cell_selection_thresh * filter_max_area
        params.t_global.cell_size = min(featureRatio, max(1, ceil(sqrt(prod(init_target_sz * search_area_scale)/(params.t_global.cell_selection_thresh * filter_max_area)))));
        featureRatio = params.t_global.cell_size;
%         featureRatio = min(featureRatio, max(1, ceil(sqrt(prod(init_target_sz * search_area_scale)/(params.t_global.cell_selection_thresh * filter_max_area)))));        
        search_area = prod(init_target_sz / featureRatio * search_area_scale);
    end
end

global_feat_params = params.t_global;

if search_area > filter_max_area
    currentScaleFactor = sqrt(search_area / filter_max_area);
else
    currentScaleFactor = 1.0;
end

% target size at the initial scale
base_target_sz = target_sz / currentScaleFactor;

% window size, taking padding into account
switch params.search_area_shape
    case 'proportional'
        sz = floor( base_target_sz * search_area_scale);     % proportional area, same aspect ratio as the target
    case 'square'
        sz = repmat(sqrt(prod(base_target_sz * search_area_scale)), 1, 2); % square area, ignores the target aspect ratio
    case 'fix_padding'
        sz = base_target_sz + sqrt(prod(base_target_sz * search_area_scale) + (base_target_sz(1) - base_target_sz(2))/4) - sum(base_target_sz)/2; % const padding
    otherwise
        error('Unknown "params.search_area_shape". Must be ''proportional'', ''square'' or ''fix_padding''');
end

% set the size to exactly match the cell size
sz = round(sz / featureRatio) * featureRatio;
use_sz = floor(sz/featureRatio);
% construct the label function- correlation output, 2D gaussian function,
% with a peak located upon the target
output_sigma = sqrt(prod(floor(base_target_sz/featureRatio))) * output_sigma_factor;
rg           = circshift(-floor((use_sz(1)-1)/2):ceil((use_sz(1)-1)/2), [0 -floor((use_sz(1)-1)/2)]);
cg           = circshift(-floor((use_sz(2)-1)/2):ceil((use_sz(2)-1)/2), [0 -floor((use_sz(2)-1)/2)]);
[rs, cs]     = ndgrid( rg,cg);
y            = exp(-0.5 * (((rs.^2 + cs.^2) / output_sigma^2)));
yf           = fft2(y); %   FFT of y.

if interpolate_response == 1
    interp_sz = use_sz * featureRatio + mod(use_sz * featureRatio-1, 2);
else
    interp_sz = use_sz + mod(use_sz-1, 2);
end

% construct cosine window
cos_window = single(hann(use_sz(1))*hann(use_sz(2))');

% Calculate feature dimension
try
    im = imread([video_path '/img/' s_frames{1}]);
catch
    try
        im = imread(s_frames{1});
    catch
        %disp([video_path '/' s_frames{1}])
        im = imread([video_path '/' s_frames{1}]);
    end
end
if size(im,3) == 3
    if all(all(im(:,:,1) == im(:,:,2)))
        colorImage = false;
    else
        colorImage = true;
    end
else
    colorImage = false;
end

% compute feature dimensionality
feature_dim = 0;
for n = 1:length(features)
    
    if ~isfield(features{n}.fparams,'useForColor')
        features{n}.fparams.useForColor = true;
    end
    
    if ~isfield(features{n}.fparams,'useForGray')
        features{n}.fparams.useForGray = true;
    end
    
    if (features{n}.fparams.useForColor && colorImage) || (features{n}.fparams.useForGray && ~colorImage)
        feature_dim = feature_dim + features{n}.fparams.nDim;
    end
end

if size(im,3) > 1 && colorImage == false
    im = im(:,:,1);
end

if nScales > 0
    scale_exp = (-floor((nScales-1)/2):ceil((nScales-1)/2));
    scaleFactors = scale_step .^ scale_exp;
    min_scale_factor = scale_step ^ ceil(log(max(5 ./ sz)) / log(scale_step));
    max_scale_factor = scale_step ^ floor(log(min([size(im,1) size(im,2)] ./ base_target_sz)) / log(scale_step));
end

if interpolate_response >= 3
    % Pre-computes the grid that is used for socre optimization
    ky = circshift(-floor((use_sz(1) - 1)/2) : ceil((use_sz(1) - 1)/2), [1, -floor((use_sz(1) - 1)/2)]);
    kx = circshift(-floor((use_sz(2) - 1)/2) : ceil((use_sz(2) - 1)/2), [1, -floor((use_sz(2) - 1)/2)])';
    newton_iterations = params.newton_iterations;
end

% initialize the projection matrix (x,y,h,w)
rect_position = zeros(num_frames, 4);
time = 0;

% allocate memory for multi-scale tracking
multires_pixel_template = zeros(sz(1), sz(2), size(im,3), nScales, 'uint8');
small_filter_sz = floor(base_target_sz/featureRatio);

% criteria = zeros(params.seq_en_frame - params.seq_st_frame,4);
loop_frame = 1;
for frame = params.seq_st_frame:params.seq_en_frame
    %load image
    try
        im = imread([video_path '/img/' s_frames{frame-params.seq_st_frame}+1]);
    catch
        try
            im = imread([s_frames{frame-params.seq_st_frame+1}]);
        catch
            im = imread([video_path '/' s_frames{frame-params.seq_st_frame+1}]);
        end
    end
    if size(im,3) > 1 && colorImage == false
        im = im(:,:,1);
    end
    if(frame == params.seq_st_frame)
        occ.total_PSR = 0;
        occ.total_Peak = 0;
        occ.avg_PSR = 0;
        occ.avg_Peak = 0;
        occ.frames = 0;
    end
    tic();
    
    %do not estimate translation and scaling on the first frame, since we
    %just want to initialize the tracker there
    if frame > params.seq_st_frame
        for scale_ind = 1:nScales
            multires_pixel_template(:,:,:,scale_ind) = ...
                get_pixels(im, pos, round(sz*currentScaleFactor*scaleFactors(scale_ind)), sz);
        end
        xtf = fft2(bsxfun(@times,get_features(multires_pixel_template,features,global_feat_params),cos_window));
        responsef = permute(sum(bsxfun(@times, conj(g_f), xtf), 3), [1 2 4 3]);
        
        % if we undersampled features, we want to interpolate the
        % response so it has the same size as the image patch
        if interpolate_response == 2
            % use dynamic interp size
            interp_sz = floor(use_sz * featureRatio * currentScaleFactor);
        end
        responsef_padded = resizeDFT2(responsef, interp_sz);
        
        % response in the spatial domain
        response = ifft2(responsef_padded, 'symmetric');
        % find maximum peak
        if interpolate_response == 3
            error('Invalid parameter value for interpolate_response');
        elseif interpolate_response == 4
            [disp_row, disp_col, sind] = resp_newton(response, responsef_padded, newton_iterations, ky, kx, use_sz);
        else
            [row, col, sind] = ind2sub(size(response), find(response == max(response(:)), 1));
            disp_row = mod(row - 1 + floor((interp_sz(1)-1)/2), interp_sz(1)) - floor((interp_sz(1)-1)/2);
            disp_col = mod(col - 1 + floor((interp_sz(2)-1)/2), interp_sz(2)) - floor((interp_sz(2)-1)/2);
        end
        %  get the potiental proposals from CF response            
        [row_prop, col_prop, val_prop, num_props] = get_proposals(response(:,:,sind),params);
        disp_row_prop = row_prop - 1 - floor((interp_sz(1)-1)/2);
        disp_col_prop = col_prop - 1 - floor((interp_sz(2)-1)/2);    
        
        % generate color likelihood_map using color_hist
        im_patch = get_pixels(im, pos, bg_area, norm_bg_area);
        % calculate color response at new predicted position
        likelihood_map = getColourMap(im_patch, bg_hist, fg_hist, params.n_bins, colorImage, cos_window_hist);
        % (TODO) in theory it should be at 0.5 (unseen colors shoud have max entropy)
        likelihood_map(isnan(likelihood_map)) = 0;    
        % each pixel of response_pwp loosely represents the likelihood that
        % the target (of size norm_target_sz) is centred on it
        [response_pwp, possible] = getCenterLikelihood(likelihood_map, norm_target_sz);
        % find the most likelihood proposals using response_pwp and
        % val_props
        [disp_row_est, disp_col_est] = proposals_estimate(disp_row_prop, disp_col_prop, val_prop, response_pwp, currentScaleFactor, area_resize_factor); 
        
        % estimate tracking status
        [Conf, occ] = tracking_status(response(:,:,sind), num_props, occ);
%         if(possible == 0)
%             possible
%         end
        % calculate translation
        switch interpolate_response
            case 0
                translation_vec = round([disp_row, disp_col] * featureRatio * currentScaleFactor * scaleFactors(sind));
                translation_vec_prop = round([disp_row_prop, disp_col_prop] * featureRatio * currentScaleFactor * scaleFactors(sind));
                translation_vec_est = round([disp_row_est, disp_col_est] * featureRatio * currentScaleFactor * scaleFactors(sind));
            case 1
                translation_vec = round([disp_row, disp_col] * currentScaleFactor * scaleFactors(sind));
                translation_vec_prop = round([disp_row_prop, disp_col_prop] * currentScaleFactor * scaleFactors(sind));
                translation_vec_est = round([disp_row_est, disp_col_est] * currentScaleFactor * scaleFactors(sind));
            case 2
                translation_vec = round([disp_row, disp_col] * scaleFactors(sind));
                translation_vec_prop = round([disp_row_prop, disp_col_prop] * scaleFactors(sind));
                translation_vec_est = round([disp_row_est, disp_col_est] * scaleFactors(sind));
            case 4
                translation_vec = round([disp_row, disp_col] * featureRatio * currentScaleFactor * scaleFactors(sind));
                translation_vec_prop = round([disp_row_prop, disp_col_prop] * featureRatio * currentScaleFactor * scaleFactors(sind));
                translation_vec_est = round([disp_row_est, disp_col_est] * featureRatio * currentScaleFactor * scaleFactors(sind));
        end
        
        % set the scale
        currentScaleFactor = currentScaleFactor * scaleFactors(sind);
        % adjust to make sure we are not to large or to small
        if currentScaleFactor < min_scale_factor
            currentScaleFactor = min_scale_factor;
        elseif currentScaleFactor > max_scale_factor
            currentScaleFactor = max_scale_factor;
        end
        
        % update position
        old_pos = pos;
        if(Conf > 2 || ~possible)
            pos = pos + translation_vec;
        else
            pos = pos + translation_vec_est;
        end
        pos_prop = old_pos + translation_vec_prop;
%         criteria(frame - params.seq_st_frame,:) = occ.criteria';
    end
        
    % extract training sample image region
    pixels = get_pixels(im,pos,round(sz*currentScaleFactor),sz);
    % patch of the target + padding
    
    % extract features and do windowing
    xf = fft2(bsxfun(@times,get_features(pixels,features,global_feat_params),cos_window));
    target_sz = floor(base_target_sz * currentScaleFactor);  
    
    if (frame == params.seq_st_frame)
        model_xf = xf;    
        model_Sxx = sum(conj(xf) .* xf, 3);
        g_f = ZA_Filter(model_xf, model_Sxx, yf, use_sz, small_filter_sz, params);
          
       % initiailize color histogram of the target
        [bg_area, fg_area, norm_bg_area, norm_target_sz, area_resize_factor] = initialHistAreas(im, target_sz, params);
        cos_window_hist = single(gausswin(norm_bg_area(1),0.5)*gausswin(norm_bg_area(2),0.5)');        
        im_patch = get_pixels(im, pos, bg_area, norm_bg_area);        
        [bg_hist, fg_hist] = initialHistModel(im_patch, bg_area, fg_area, target_sz, norm_bg_area, params.n_bins, colorImage);
    else
        if(Conf > 1.93)
            model_xf = ((1 - learning_rate) * model_xf) + (learning_rate * xf);       
            model_Sxx = ((1 - learning_rate) * model_Sxx) + (learning_rate * sum(conj(xf) .* xf, 3));
            g_f = ZA_Filter(model_xf, model_Sxx, yf, use_sz, small_filter_sz, params);      

            % update color histogram of the target
            [bg_area, fg_area, area_resize_factor] = updateHistAreas(im, target_sz, norm_bg_area, params);
            im_patch = get_pixels(im, pos, bg_area, norm_bg_area); 
            [bg_hist, fg_hist] = updateHistModel(im_patch, bg_area, fg_area, target_sz, norm_bg_area, bg_hist, fg_hist, params.n_bins, colorImage, params.learning_rate_pwp);
        end
    end    
    
    %save position and calculate FPS
    rect_position(loop_frame,:) = [pos([2,1]) - floor(target_sz([2,1])/2), target_sz([2,1])];
    time = time + toc();
    
    %visualization
    if visualization == 1
        rect_position_vis = [pos([2,1]) - target_sz([2,1])/2, target_sz([2,1])];
        im_to_show = double(im)/255;
        if size(im_to_show,3) == 1
            im_to_show = repmat(im_to_show, [1 1 3]);
        end
        if frame == params.seq_st_frame
            fig_handle = figure('Name', 'Tracking');
%             q=get(fig_handle,'position');
%             q(1)=0;%ÉèÖÃ×ó±ß¾àÀëÖµÎªÁã
%             q(2)=0;%ÉèÖÃÓÒ±ß¾àÀëÖµÎªÁã
%             set(fig_handle, 'Units', 'pixels', 'Position', q);         
            imagesc(im_to_show);
            hold on;
            rectangle('Position',rect_position_vis, 'EdgeColor','g', 'LineWidth',2);
            text(10, 10, int2str(frame), 'color', [0 1 1]);
            axis image;set(gca, 'Units', 'normalized', 'Position', [0 0 1 1]);
            hold off;
        else

%             resp_sz = round(sz*currentScaleFactor*scaleFactors(scale_ind));
%             xs = floor(old_pos(2)) + (1:resp_sz(2)) - floor(resp_sz(2)/2);
%             ys = floor(old_pos(1)) + (1:resp_sz(1)) - floor(resp_sz(1)/2);
%             sc_ind = floor((nScales - 1)/2) + 1;
%             rect_positon_est = [pos_est([2,1]) - target_sz([2,1])/2, target_sz([2,1])];
            [N,~] = size(pos_prop);
            rect_position_prop = [pos_prop(:,[2,1]) - target_sz([2,1])/2, repmat(target_sz([2,1]),N,1)];          
            figure(fig_handle);
            imagesc(im_to_show);
            hold on;
%             resp_handle = imagesc(xs, ys, fftshift(response(:,:,sc_ind))); colormap hsv;
%             alpha(resp_handle, 0.2);         
            for i = 1:N 
                rectangle('Position',rect_position_prop(i,:), 'EdgeColor','b', 'LineWidth',2);
            end
            rectangle('Position',rect_position_vis, 'EdgeColor','g', 'LineWidth',2); 
%             rectangle('Position',rect_positon_est, 'EdgeColor','r', 'LineWidth',2);
            text(20, 20, ['# Frame : ' int2str(frame) ' / ' int2str(params.seq_en_frame)], 'color', [1 0 0], 'BackgroundColor', [1 1 1], 'fontsize', 16);
            text(20, 40, ['FPS : ' num2str(1/(time/loop_frame))], 'color', [1 0 0], 'BackgroundColor', [1 1 1], 'fontsize', 16);
            text(10, 10, int2str(frame), 'color', [0 1 1]);
            axis image;set(gca, 'Units', 'normalized', 'Position', [0 0 1 1]);             
            hold off;
        end
        drawnow
%         screen = getframe(fig_handle,[0,0,560,420]);
%         im_save = frame2im(screen);
%         filename = [video_path, '/results/', num2str(frame,'%04d'), '.jpg']; 
%         imwrite(im_save,filename,'jpg');
    end
    if debug == 1 && frame > params.seq_st_frame
        im_patch = get_pixels(im, old_pos, bg_area, norm_bg_area);
        pixels = get_pixels(im,old_pos,round(sz*currentScaleFactor),sz);      
        mySubplot(2,2,3,1,im_patch,'FG+BG','gray');
        mySubplot(2,2,3,2,likelihood_map,'obj.likelihood','parula'); 
        mySubplot(2,2,3,3,response_pwp,'obj.likelihood','parula'); 
        mySubplot(2,2,3,4,pixels,'Search area','gray');
        mySubplot(2,2,3,5,fftshift(response(:,:,sind)),'CF likelihood','parula');
    end
 
    loop_frame = loop_frame + 1;
    
end
%   save resutls.
fps = loop_frame / time;
results.type = 'rect';
results.res = rect_position;
results.fps = fps;
% results.criteria = criteria;
