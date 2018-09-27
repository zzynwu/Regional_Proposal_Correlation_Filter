
clc;clear;
obj = VideoReader('D:/Phd/dataset/valid_BW.avi');%输入视频位置
numFrames = obj.NumberOfFrames;% 帧的总数

seq.len = numFrames;
seq.init_rect = [805 445 55 21];%[805 445 56 23]  %[873 13 47 19]  %[993 315 50 18]
seq.VidName = 'Car';
seq.st_frame = 5057;  %5057 %3662 %100
seq.en_frame = seq.len;  %seq.len %5056 %3521
seq.last_frame = seq.len;

img_path = 'D:/phd/dataset/Carwhite';
img_files = dir(fullfile(img_path, '/img/*.jpg'));
img_files = {img_files.name};
seq.s_frames = cellstr(img_files);

%   Run BACF- main function
learning_rate = 0.013;  %   you can use different learning rate for different benchmarks.
results       = run_Tracker(seq, img_path, learning_rate);