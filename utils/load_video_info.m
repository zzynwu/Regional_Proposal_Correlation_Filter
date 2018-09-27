function [seq, ground_truth] = load_video_info(video_path)

if exist([video_path '/groundtruth_rect.txt'],'file')
    ground_truth = dlmread([video_path '/groundtruth_rect.txt']);
elseif exist([video_path '/groundtruth_rect.1.txt'],'file')
    ground_truth = dlmread([video_path '/groundtruth_rect.1.txt']);
elseif exist([video_path '/groundtruth_rect.2.txt'],'file')
    ground_truth = dlmread([video_path '/groundtruth_rect.2.txt']);
end

seq.len = size(ground_truth, 1);
seq.init_rect = ground_truth(1,:);% [270,311,22,65]

img_path = [video_path '/img/'];

img_files = dir(fullfile(img_path, '*.jpg'));
% img_files = img_files(4:end);
img_files = {img_files.name};
% img_files = [img_path img_files];
% if exist([img_path num2str(1, '%04i.png')], 'file'),
%     img_files = num2str((1:seq.len)', [img_path '%04i.png']);
% elseif exist([img_path num2str(1, '%04i.jpg')], 'file'),
%     img_files = num2str((1:seq.len)', [img_path '%04i.jpg']);
% elseif exist([img_path num2str(1, '%04i.bmp')], 'file'),
%     img_files = num2str((1:seq.len)', [img_path '%04i.bmp']);
% else
%     error('No image files to load.')
% end

seq.s_frames = cellstr(img_files);

end

