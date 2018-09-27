clc
clear
framesPath = 'D:/Phd/dataset/Carwhite/results';%ͼ����������·����ͬʱҪ��֤ͼ���С��ͬ
videoName = 'Carwhite.avi';%��ʾ��Ҫ��������Ƶ�ļ�������
fps = 25; %֡��

img_files = dir(fullfile(framesPath, '*.jpg'));
% img_files = img_files(4:end);
img_files = {img_files.name};
n = length(img_files);
startFrame = 1; %����һ֡��ʼ
endFrame = n; %��һ֡����

if(exist('videoName','file'))
    delete videoName.avi
end

%������Ƶ�Ĳ����趨
aviobj=VideoWriter(videoName);  %����һ��avi��Ƶ�ļ����󣬿�ʼʱ��Ϊ��
aviobj.FrameRate=fps;

open(aviobj);%Open file for writing video data
%����ͼƬ
for i=startFrame:endFrame
    fileName = img_files{1,1};
    frames=imread([framesPath,'/' ,fileName]);
    writeVideo(aviobj,frames);
end
close(aviobj);% �رմ�����Ƶ