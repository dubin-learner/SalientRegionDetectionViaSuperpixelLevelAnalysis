clear; close all; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Du Bin in Beijing Institute of Technology (BIT)
% Last modified time : 2016/5/22
% Paper: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imgRoot = './InputImage/';		% Test Image Path
supDir  = './Superpixel/';      % Superpixel Segmentation Results

addpath(imgRoot);
addpath('Functions');
mkdir(supDir);
mkdir('./DownSample/');
mkdir('./Results/');

imageNum    = 1;            % Number of Test Images
downRate    = 4;            % Downsampling Rate
spNum       = 400;          % Number of Superpixel you want
backsRatio  = 0.3;          % Background Proportions

for img_i = 1:imageNum
%% Input Image and Downsample
tic;
imName      = [num2str(img_i), '.bmp'];
img_raw     = imread(imName);
img_down    = img_raw(1:downRate:end, 1:downRate:end, :);

imwrite(img_down, ['./DownSample/', imName]);

img_down    = double(img_down)/255;
img_down    = colorspace('Lab<-RGB', img_down);
img_gray    = img_down(:, :, 1);

[row, col, ~] = size(img_down);

%% Generate Superpixels
supDir  = './Superpixel/';
cmd     = ['SLICSuperpixelSegmentation', ' ', ...
    ['./DownSample/', imName], ' ', '20', ' ', int2str(spNum), ' ', supDir];
system(cmd);

superpixels = ReadDAT([row, col], [supDir, imName(1:end-4), '.dat']);
spNum       = max( superpixels(:) );    % Number of Superpixel in fact

adjloop     = AdjcProcloop(superpixels, spNum);
% toc;

%% The representation value of a superpixel
% The Centroid of every superpixels
S = regionprops(superpixels, 'Centroid');
positions = zeros(spNum, 2);

% value of every superpixels
mean_value = zeros(1, spNum);
for sp_i = 1:spNum
    region_i = superpixels == sp_i;
    mean_value(sp_i)    = mean( img_gray(region_i) );
    positions(sp_i, :)  = S(sp_i).Centroid;
end
mean_value = mean_value/max(mean_value);

%% Gradient of the representation value
grad_x = zeros(1, spNum);
grad_y = zeros(1, spNum);
grad_g = zeros(1, spNum);

sigma = 1;
for sp_i = 1:spNum
    neigbors = find( adjloop(sp_i, :) == 1 );
    % Centroid Distance between every neigbor and current superpixel
    elements = positions(neigbors, :) - ones(length(neigbors), 1)*positions(sp_i, :);
    % Gradient divided by the distance
    grad_sum = ( mean_value(neigbors) - mean_value(sp_i) ).*...
        exp(- ( sqrt( elements(:, 1).^2 + elements(:, 2).^2 ) )/sigma )';
    
    elements_x = elements(:, 1)'.*grad_sum;
    elements_y = elements(:, 2)'.*grad_sum;
    grad_x(sp_i) = sum(elements_x);
    grad_y(sp_i) = sum(elements_y);
    
    tmp4 = sqrt(elements(:, 1).^2 + elements(:, 2).^2);   
    grad_g(sp_i) = mean_value(sp_i) + ...
        mean( 1/sqrt(2*pi)*exp( -tmp4.^2/2/sigma^2 ) );
end

grad_x = grad_x/max(grad_x);
grad_y = grad_y/max(grad_y);
grad_g = grad_g/max(grad_g);

%% Structure Tensor for Texture Saliency Map
textureSal_sp = zeros(1, spNum);
for sp_i = 1:spNum
    neigbors    = find( adjloop(sp_i, :) == 1 );
    vectors     = [neigbors, sp_i];
    
    matrix = [  sum(grad_x(vectors).^2), sum(grad_x(vectors).*grad_y(vectors)), sum(grad_x(vectors).*grad_g(vectors));
                sum(grad_y(vectors).*grad_x(vectors)), sum(grad_y(vectors).^2), sum(grad_y(vectors).*grad_g(vectors));
                sum(grad_g(vectors).*grad_x(vectors)), sum(grad_g(vectors).*grad_y(vectors)), sum(grad_g(vectors).^2)];
%     matrix = [  grad_x(sp_i).^2, grad_x(sp_i)*grad_y(sp_i), grad_x(sp_i)*grad_g(sp_i);
%                 grad_x(sp_i)*grad_y(sp_i), grad_y(sp_i).^2, grad_y(sp_i)*grad_g(sp_i);
%                 grad_x(sp_i)*grad_g(sp_i), grad_y(sp_i)*grad_g(sp_i), grad_g(sp_i).^2];

    matrix = matrix/length(vectors);
    
    [~, tez] = eig(matrix);
    tez = diag(tez);
    
    textureSal_sp(sp_i) = ( max(tez) - min(tez) )^2/sum(tez);

end
textureSal_sp = ( textureSal_sp - min(textureSal_sp) )/...
    ( max(textureSal_sp) - min(textureSal_sp) );

%% Calculate Color Distance for Color Saliency Map
% Color Center of each superpixels
cls = zeros(1, spNum);
for sp_i = 1:spNum
    cls(sp_i) = sp_i;
end
colorCenter = computeColorCenter(img_down, superpixels, cls, spNum);

backSpNum   = round(spNum*backsRatio);
spSalSorted = sort(textureSal_sp);
backSp      = find( textureSal_sp <= spSalSorted(backSpNum) );
backSpValue = mean( colorCenter(backSp, :), 1 );
colorDist   = colorCenter - ones(spNum, 1)*backSpValue;
colorSal_sp = sum(colorDist.^2, 2);
colorSal_sp = colorSal_sp'/max(colorSal_sp);

%% Fusion and re-scale
salMap_s    = textureSal_sp.*colorSal_sp;
spResize    = imresize(superpixels, downRate, 'nearest');

img_raw     = imread(imName);
img_raw     = rgb2gray(img_raw);

% toc;
%% Superpixel level to Pixel level
salMap = sp2pixel(img_raw, spNum, spResize, adjloop, salMap_s);
supDir = './Results/';
imwrite(salMap, [supDir, imName(1:end-4), '.tif']);

toc;
%% Display Result
figure, imshow(salMap, []);
title(['Salient Region Detection Result ', num2str(img_i)]);

end

