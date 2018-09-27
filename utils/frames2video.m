clc
clear
framesPath = 'D:/Phd/dataset/Carwhite/results';%图像序列所在路径，同时要保证图像大小相同
videoName = 'Carwhite.avi';%表示将要创建的视频文件的名字
fps = 25; %帧率

img_files = dir(fullfile(framesPath, '*.jpg'));
% img_files = img_files(4:end);
img_files = {img_files.name};
n = length(img_files);
startFrame = 1; %从哪一帧开始
endFrame = n; %哪一帧结束

if(exist('videoName','file'))
    delete videoName.avi
end

%生成视频的参数设定
aviobj=VideoWriter(videoName);  %创建一个avi视频文件对象，开始时其为空
aviobj.FrameRate=fps;

open(aviobj);%Open file for writing video data
%读入图片
for i=startFrame:endFrame
    fileName = img_files{1,1};
    frames=imread([framesPath,'/' ,fileName]);
    writeVideo(aviobj,frames);
end
close(aviobj);% 关闭创建视频