function [DICpara] = ReadImageMaskHole(DICpara,file_name)
%FUNCTION [DICpara] = ReadImageMask(DICpara,file_name)
% ----------------------------------------------
%   This script is to generate a binary image mask in the reference image
%   This code only works for a single hole in a 2D sample
%   For other general geometries, please see code "ReadImageMask.m"
%  
%   INPUT: DICpara       DIC parameters
%          file_name     Loaded DIC raw images
%
%   OUTPUT: DICpara      DIC parameters
%
% ----------------------------------------------
% Author: Jin Yang.  
% Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
% Last time updated: 12/2020.
% ==============================================

 
%%
% ====== Generate an image mask in the reference image ======
ImgRef = (imread(file_name{1})); 
try ImgRef = rgb2gray(ImgRef); catch; end % Change the image to grayscale image if it is rgb image

%% ====== Find center of the hole ======
figure, imshow(ImgRef);  % Read reference image
title('Click more than 5 points along the hole boundary, then press -Enter- key','fontweight','normal');
fprintf('\n'); fprintf('--- Define hole geometry ---  \n');
fprintf('Click more than 5 points along the hole boundary,  \n');
fprintf('then press -Enter- key. \n');
[gridx,gridy] = ginput(); ParHole = CircleFitByTaubin([gridx,gridy]); % ParHole is [xCenCoord, yCenCoord, rad]; ParHole=[194.146,542.51,52.072];
fprintf('--- Define hole geometry: Done! ---  \n');

 
%% ====== Construct the image mask ======
ImgRefMask = ones(size(ImgRef));
[tempxx,tempyy] = ndgrid(1:size(ImgRefMask,1), 1:size(ImgRefMask,2));
dist2HoleCenter = sqrt( (tempxx-ParHole(2)).^2 + (tempyy-ParHole(1)).^2 ); % Compute the distance between image points to the hole center

[row,col] = find(dist2HoleCenter < ParHole(3));
for tempi = 1:length(row)
   ImgRefMask(row(tempi),col(tempi)) = 0; 
end

% Plot
figure, imshow(ImgRefMask); title('Image Mask')

removeobjradius =  round(0.2*ParHole(3)); % remove all object containing fewer than hole Rad
se = strel('disk',removeobjradius); ImgRefMask = imclose(ImgRefMask,se); imshow(ImgRefMask);


%% Deal with outside  DICpara.gridxyROIRange
gridx = DICpara.gridxyROIRange.gridx;
gridy = DICpara.gridxyROIRange.gridy;

ImgRefMask(:, 1:gridx(1)) = 1;
ImgRefMask(:, gridx(2):end) = 1;
ImgRefMask(1:gridy(1), :) = 1;
ImgRefMask(gridy(2):end, :) = 1;

% ===== Plot =====
close all; figure,
imshow(ImgRefMask); title('Ref image mask','interpreter','latex'); 
axis equal; axis tight; set(gca,'fontsize',18); set(gcf,'color','w'); box on;
xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
a = gca; a.TickLabelInterpreter = 'latex'; 


%% ====== Store ImgRefMask ======
DICpara.ImgRefMask = ImgRefMask'; % Consider the image coordinates
DICpara.ParHole = ParHole;



end