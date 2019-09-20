clc; clear all; close all;

addpath(genpath("./YUV2Image"));

fnames = [
    "akiyo_qcif.yuv";
    "carphone_qcif.yuv";
    "container_qcif.yuv";
    "foreman_qcif.yuv";
    "suzie_qcif.yuv";
];

width = 176;
height = 144;

X = [];
y = [];

for i = 1:length(fnames)
    
    fileName = fnames(i);
    
    fprintf('Loading video %s... ', fileName)
    
    idxFrame = 1;
    while idxFrame < 1000
        try
            [mov,imgRgb] = loadFileYuv(fileName, width, height, idxFrame);
            
            imgDouble = im2double(imgRgb);
            imgGrey = rgb2gray(imgDouble);
            imgFlat = reshape(imgGrey, 1, width*height);
            
            X = [X; imgFlat];
            y = [y; i];
            
            idxFrame = idxFrame + 1;
        catch ME
            fprintf('finished!\n')
            break
        end
    end
    
end

save('./X_video.mat','X');
save('./y_video.mat','y');
