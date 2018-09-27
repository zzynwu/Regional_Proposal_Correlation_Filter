clc
clear
obj = VideoReader('D:\Phd\dataset\valid_BW.avi');%输入视频位置
numFrames = obj.NumberOfFrames;% 帧的总数
 for k = 1 : numFrames% 读取前15帧
     frame = read(obj,k);%读取第几帧
    % imshow(frame);%显示帧
     filename = ['D:\phd\dataset\Carwhite\img\',num2str(k,'%04d'),'.jpg'];
     imwrite(frame,filename,'jpg');% 保存帧
 end