clc
clear
obj = VideoReader('D:\Phd\dataset\valid_BW.avi');%������Ƶλ��
numFrames = obj.NumberOfFrames;% ֡������
 for k = 1 : numFrames% ��ȡǰ15֡
     frame = read(obj,k);%��ȡ�ڼ�֡
    % imshow(frame);%��ʾ֡
     filename = ['D:\phd\dataset\Carwhite\img\',num2str(k,'%04d'),'.jpg'];
     imwrite(frame,filename,'jpg');% ����֡
 end