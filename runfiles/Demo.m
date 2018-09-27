
%   This script runs the original implementation of Background Aware Correlation Filters (BACF) for visual tracking.
%   the code is tested for Mac, Windows and Linux- you may need to compile
%   some of the mex files.
%   Paper is published in ICCV 2017- 
%   Some functions are borrowed from other papers (SRDCF, CCOT, KCF, etc)- and
%   their copyright belongs to the paper's authors.
%   copyright- Hamed Kiani
%   contact me: hamedkg@gmail.com


%   This demo runs on OTB50, you can use any benchmark by setting the seq
%   path, and using the standard annotation txt files.
clear;clc;
close all;
% Load video information
% base_path  = '../BACF_ICCV_paper/sequences/Benchmark_TB_100';
base_path = 'D:/Phd/dataset/Benchmark2015';
% videos = {'Basketball', 'Bolt', 'Boy', 'Car4', 'CarDark', 'CarScale', ...
%     'Coke', 'Couple', 'Crossing', 'David2', 'David3', 'David', 'Deer', ...
%     'Dog1', 'Doll', 'Dudek', 'Faceocc1', 'Faceocc2', 'Fish', 'Fleetface', ...
%     'Football', 'Football1', 'Freeman1', 'Freeman3', 'Freeman4', 'Girl', ...
%     'Ironman', 'Jogging', 'Jumping', 'Lemming', 'Liquor', 'Matrix', ...
%     'Mhyang', 'MotorRolling', 'MountainBike', 'Shaking', 'Singer1', ...
%     'Singer2', 'Skating1', 'Skiing', 'Soccer', 'Subway', 'Suv', 'Sylvester', ...
%     'Tiger1', 'Tiger2', 'Trellis', 'Walking', 'Walking2', 'Woman',};
videos = {'Basketball', 'Bolt', 'Boy', 'Car4', 'CarDark', 'CarScale', ...
    'Coke', 'Couple', 'Crossing', 'David2', 'David3', 'David', 'Deer', ...
    'Dog1', 'Doll', 'Dudek', 'Faceocc1', 'Faceocc2', 'Fish', 'Fleetface', ...
    'Football', 'Football1', 'Freeman1', 'Freeman3', 'Freeman4', 'Girl', ...
    'Ironman', 'Jogging', 'Jumping', 'Lemming', 'Liquor', 'Matrix', ...
    'Mhyang', 'MotorRolling', 'MountainBike', 'Shaking', 'Singer1', ...
    'Singer2', 'Skating1', 'Skiing', 'Soccer', 'Subway', 'Suv', 'Sylvester', ...
    'Tiger1', 'Tiger2', 'Trellis', 'Walking', 'Walking2', 'Woman',...
    'Biker','Bird1','Bird2','Blurbody','BlurCar1','BlurCar2',...
    'BlurCar3','BlurCar4','BlurFace','BlurOwl','Board','box',...
    'bolt2','car1','car2','car24','Clifbar','Coupon','Crowds',...
    'Dancer','Dancer2','Diving','Dog','Dragonbaby','Girl2',...
    'Gym','human2','human3','human4','human5','human6','human7',...
    'human8','human9','Jump','KiteSurf','Man','Panda','RedTeam',...
    'Rubik','skater','skater2','skating2','Surfer','Toy','Trans',...
    'Twinnings','Vase',};
OPs = zeros(numel(videos),1);
FPSs = zeros(numel(videos),1);
OPs_OTB50 = [];
FPSs_OTB50 = [];

for vid = 1:numel(videos)
    close all;
    video_path = [base_path '/' videos{vid}];
    [seq, ground_truth] = load_video_info(video_path);
    seq.VidName = videos{vid};
%     if(exist([video_path '/criteria.txt'],'file'))
%         continue;
%     end
    st_frame = 1;
    en_frame = seq.len;
    if (strcmp(videos{vid}, 'David'))
        st_frame = 300;
        en_frame = 770;
    elseif (strcmp(videos{vid}, 'Football1'))
        st_frame = 1;
        en_frame = 74;
    elseif (strcmp(videos{vid}, 'Freeman3'))
        st_frame = 1;
        en_frame = 460;
    elseif (strcmp(videos{vid}, 'Freeman4'))
        st_frame = 1;
        en_frame = 283;
    end
    seq.st_frame = st_frame;
    seq.en_frame = en_frame;
    gt_boxes = [ground_truth(:,1:2), ground_truth(:,1:2) + ground_truth(:,3:4) - ones(size(ground_truth,1), 2)];
    
    % Run BACF- main function
    learning_rate = 0.013;  %   you can use different learning rate for different benchmarks.
    results = run_Tracker(seq, video_path, learning_rate);
%     fid = fopen([video_path '/criteria_interp.txt'] ,'w');
%     fprintf( fid, '%f \r\n', results.criteria);
%     fclose(fid);
    results.gt = gt_boxes;
    %   compute the OP
    pd_boxes = results.res;
    pd_boxes = [pd_boxes(:,1:2), pd_boxes(:,1:2) + pd_boxes(:,3:4) - ones(size(pd_boxes,1), 2)  ];
    OP = zeros(size(gt_boxes,1),1);
    for i=1:size(gt_boxes,1)
        b_gt = gt_boxes(i,:);
        b_pd = pd_boxes(i,:);
        OP(i) = computePascalScore(b_gt,b_pd);
    end
    OPs(vid) = sum(OP >= 0.5) / numel(OP);
    FPSs(vid) = results.fps;
    display([videos{vid}  '---->' '   FPS:   ' num2str(results.fps)   '    op:   '   num2str(OPs(vid))]);
   
    FPSs_OTB50 = [FPSs_OTB50; results.fps];
    OPs_OTB50 =  [OPs_OTB50; OPs(vid)];
end
