function [bg_area, fg_area, area_resize_factor] = updateHistAreas(im, cur_target_sz, norm_bg_area, params)	
    % we want a regular frame surrounding the object
	avg_dim = sum(cur_target_sz); %sum(cur_target_sz)/2;
	% size from which we extract features
	bg_area = round(cur_target_sz + avg_dim);
	% pick a "safe" region smaller than bbox to avoid mislabeling
	fg_area = round(cur_target_sz - avg_dim * params.inner_padding);
	% saturate to image size
	if(bg_area(2)>size(im,2)), bg_area(2)=size(im,2)-1; end
	if(bg_area(1)>size(im,1)), bg_area(1)=size(im,1)-1; end
	% make sure the differences are a multiple of 2 (makes things easier later in color histograms)
	bg_area = bg_area - mod(bg_area - cur_target_sz, 2);
	fg_area = fg_area + mod(bg_area - fg_area, 2);
    
    % Compute the rectangle with (or close to) params.fixedArea and
    % same aspect ratio as the target bbox
	area_resize_factor = sqrt(prod(norm_bg_area)/prod(bg_area));
% 	norm_bg_area = round(bg_area * area_resize_factor);
%     norm_target_sz = round(cur_target_sz * area_resize_factor);
end
