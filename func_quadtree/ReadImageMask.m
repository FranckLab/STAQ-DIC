function [DICpara] = ReadImageMask(DICpara)
%FUNCTION [DICpara] = ReadImageMask(DICpara)
% ----------------------------------------------
%   This script is to generate a binary image mask in the reference image
%  
%   INPUT: DICpara       DIC parameters
%
%   OUTPUT: DICpara      DIC parameters
%
% ----------------------------------------------
% Author: Jin Yang.  
% Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
% Last time updated: 12/2020.
% ==============================================
 
 
%%
% ====== Generate an image mask for the reference image ======
disp('--- Please load a reference frame to generate an image binary mask ---')
file_name_ImgMask{1,1} = uigetfile('*.tif','Select the reference frame to generate an image mask');
        
ImgRef = (imread(file_name_ImgMask{1,1})); 
try ImgRef = rgb2gray(ImgRef(:,:,1:3)); catch; end % Change the image to grayscale image if it is rgb image
disp('Reference image is loaded.')

 
%% ====== Get the grayscale value threshold and construct the image mask ======
ImgRefGaussFilt = imgaussfilt(imgaussfilt(imgaussfilt(ImgRef,1),1),1); % Gaussian filter

close all; figure, imshow(ImgRefGaussFilt);  % Read reference image
title('Click more than 6 points in the background, then press -Enter-','fontweight','normal');
fprintf('\n');
fprintf('--- Define image mask file step 1: background ---  \n');
fprintf('Click more than 6 points in the background, \n');
fprintf('then press the -Enter- key. \n');
[gridx,gridy] = ginput();  ImgRefMaskThresholdInside = 0;
for tempi = 1:length(gridx)
    ImgRefMaskThresholdInside = ImgRefMaskThresholdInside + double(ImgRefGaussFilt(round(gridy(tempi)),round(gridx(tempi))));
end
ImgRefMaskThresholdInside = ImgRefMaskThresholdInside/length(gridx);
% For example, ImgRefMaskThreshold = 20; in the imgFiles_Quadtree_Sample12 exampple

close all; figure, imshow(ImgRefGaussFilt);  % Read reference image
title('Click more than 6 points on the sample pattern, then press -Enter-','fontweight','normal');
fprintf('\n');
fprintf('--- Define image mask file step 2: sample surface ---  \n');
fprintf('Click more than 6 points on the sample pattern, \n');
fprintf('then press the -Enter- key. \n');
[gridx,gridy] = ginput();  ImgRefMaskThresholdOutside = 0;
for tempi = 1:length(gridx)
    ImgRefMaskThresholdOutside = ImgRefMaskThresholdOutside + double(ImgRefGaussFilt(round(gridy(tempi)),round(gridx(tempi))));
end
ImgRefMaskThresholdOutside = ImgRefMaskThresholdOutside/length(gridx);
% For example, ImgRefMaskThreshold = 20; in the imgFiles_Quadtree_Sample12 exampple

fprintf('\n'); fprintf('Is the background darker or brighter than the sample surface?  \n');
fprintf('    0: Background is darker; \n');
fprintf('    1: Background is brighter; \n');
prompt = 'Input here: '; HoleDarkerOrBrighter = input(prompt);

if HoleDarkerOrBrighter == 0
    ImgRefMask = logical(ImgRefGaussFilt > 0.5*(ImgRefMaskThresholdInside+ImgRefMaskThresholdOutside));
elseif HoleDarkerOrBrighter == 1
    ImgRefMask = logical(ImgRefGaussFilt < 0.5*(ImgRefMaskThresholdInside+ImgRefMaskThresholdOutside));
end

% [tempxx,tempyy] = ndgrid(1:size(ImgRefMask,1), 1:size(ImgRefMask,2));
% dist2HoleCenter = sqrt( (tempxx-ParHole(2)).^2 + (tempyy-ParHole(1)).^2 );
% [row,col] = find(dist2HoleCenter < ParHole(3));
% for tempi = 1:length(row)
%    ImgRefMask(row(tempi),col(tempi)) = 0; 
% end
figure, imshow(ImgRefMask); title('Image binary mask')


RemoveObjRadiusOKOrNot = 1; ImgRefMasktemp=ImgRefMask;
while RemoveObjRadiusOKOrNot == 1
    fprintf('\n');
    fprintf('Remove all objects smaller than the critical radius. \n');
    prompt = 'Input here: ';
    removeobjradius = round(input(prompt));
    se = strel('disk',removeobjradius); ImgRefMasktemp = imclose(ImgRefMask,se); 
    close all; figure, imshow(ImgRefMasktemp);
    fprintf('Are you satisfied with current removal radius? (0-Yes; 1-No) \n');
    prompt = 'Input here: ';
    RemoveObjRadiusOKOrNot = input(prompt);
end



% removeobjradius =  round(0.2*ParHole(3)); % remove all object containing fewer than hole Rad
se = strel('disk',removeobjradius); ImgRefMask = imclose(ImgRefMask,se); imshow(ImgRefMask);


%% Deal with outside  DICpara.gridxyROIRange
gridx = DICpara.gridxyROIRange.gridx;
gridy = DICpara.gridxyROIRange.gridy;

ImgRefMask(:, 1:gridx(1)) = 1;
ImgRefMask(:, gridx(2):end) = 1;
ImgRefMask(1:gridy(1), :) = 1;
ImgRefMask(gridy(2):end, :) = 1;

imshow(ImgRefMask);
 

%% ====== Store ImgRefMask ======
DICpara.ImgRefMask = ImgRefMask'; % Consider the image coordinates


end
