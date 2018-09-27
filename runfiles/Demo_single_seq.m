
%   This script runs the original implementation of Background Aware Correlation Filters (BACF) for visual tracking.
%   the code is tested for Mac, Windows and Linux- you may need to compile
%   some of the mex files.
%   Paper is published in ICCV 2017- Italy
%   Some functions are borrowed from other papers (SRDCF, CCOT, KCF, etc)- and
%   their copyright belongs to the paper's authors.
%   copyright- Hamed Kiani (CMU, RI, 2017)

%   contact me: hamedkg@gmail.com

clear; clc; close all;
%   Load video information
base_path = 'D:/Phd/dataset/Benchmark2015';
video      = 'girl';

video_path = [base_path '/' video];
[seq, ground_truth] = load_video_info(video_path);
seq.VidName = video;
seq.st_frame = 1;
seq.en_frame = seq.len;
    if (strcmp(video, 'David'))
        st_frame = 300;
        en_frame = 770;
    elseif (strcmp(video, 'Football1'))
        st_frame = 1;
        en_frame = 74;
    elseif (strcmp(video, 'Freeman3'))
        st_frame = 1;
        en_frame = 460;
    elseif (strcmp(video, 'Freeman4'))
        st_frame = 1;
        en_frame = 283;
    end
gt_boxes = [ground_truth(:,1:2), ground_truth(:,1:2) + ground_truth(:,3:4) - ones(size(ground_truth,1), 2)];

%   Run BACF- main function
learning_rate = 0.013;  %   you can use different learning rate for different benchmarks.
results       = run_Tracker(seq, video_path, learning_rate);
fid = fopen([video_path '/criteria.txt'] ,'w');
fprintf( fid, '%f %f %f %d \r\n', results.criteria');
fclose(fid);
%   compute the OP
pd_boxes = results.res;
pd_boxes = [pd_boxes(:,1:2), pd_boxes(:,1:2) + pd_boxes(:,3:4) - ones(size(pd_boxes,1), 2)  ];
OP = zeros(size(gt_boxes,1),1);
for i=1:size(gt_boxes,1)
    b_gt = gt_boxes(i,:);
    b_pd = pd_boxes(i,:);
    OP(i) = computePascalScore(b_gt,b_pd);
end
OP_vid = sum(OP >= 0.5) / numel(OP);
FPS_vid = results.fps;
display([video  '---->' '   FPS:   ' num2str(FPS_vid)   '    op:   '   num2str(OP_vid)]);
