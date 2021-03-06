clear; clc;
close all;

addpath utils;
addpath 'qcut_utils';

% parent_path = 'simupor_dataset/';
% load([parent_path 'datasets.mat']);
% total_columns = 1;

inputBrain = '/rsrch1/ip/rglenn1/quantumSegmentation/QIS-Net/LiverInput/13.000000-t1vibeqfstrap2bhFIL-72776';
outputFolder = '/rsrch1/ip/rglenn1/quantumSegmentation/qcuts3d/LiverOutput/liver';
dicomlist = dir(fullfile(inputBrain,'*.dcm'));
total_Im = numel(dicomlist);
disp(total_Im); 
%[V, spatial, dim] = dicomreadVolume(fullfile(matlabroot,inputBrain));
filename = fullfile(inputBrain,dicomlist(1).name);

img = dicomread(filename);
imgSize = size(img);
width = imgSize(1);
height = imgSize(2);

dcm3d = zeros( width, height, total_Im, 'uint8');

for cnt = 1 : total_Im 
    disp('cnt');
    disp(cnt);
    filename = fullfile(inputBrain,dicomlist(cnt).name);
    img = dicomread(filename);
 
    % Normalize data
    img = double(img)/255;
    dcm3d(:,:,cnt) = img;
end

    
% Apply QCUT
SalMap = applyQCUTv3(dcm3d);

output_mean=mean(SalMap,4);
output = double(output_mean>=0.75);
for cnt = 1: total_Im
    ouput = output_mean(:,:,cnt);
    filename = sprintf('LiverOutput/liver_%d.png', cnt);
    imshow(ouput);
    biggerIm = resizeImage(ouput* 255, 5/6);
    imwrite(biggerIm, filename);
end

function biggerIm = resizeImage(img, scaleFactor)
    sz = size(img);
    xg = 1:sz(1);
    yg = 1:sz(2);
    F = griddedInterpolant({xg,yg},double(img));
    xq = (0:scaleFactor:sz(1))';
    yq = (0:scaleFactor:sz(2))';
    biggerIm = uint8(F({xq, yq}));
   
end
  