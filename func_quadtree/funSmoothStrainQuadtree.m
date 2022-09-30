function F = funSmoothStrainQuadtree(F,DICmesh,DICpara)
%FUNSMOOTHSTRAINQUADTREE: to smooth solved strain fields by curvature regularization
% 	F = funSmoothStrainQuadtree(F,DICmesh,DICpara)
% ----------------------------------------------
%
%   INPUT: F                 Deformation gradient tensor: 
%                            F = [F11_node1, F21_node1, F12_node1, F22_node1, ... , F11_nodeN, F21_nodeN, F12_nodeN, F22_nodeN]';
%          DICmesh           DIC mesh
%          DICpara           DIC parameters
%
%   OUTPUT: F                Smoothed strain fields by curvature regularization
%
% ----------------------------------------------
% Reference
% [1] RegularizeNd. Matlab File Exchange open source. 
% https://www.mathworks.com/matlabcentral/fileexchange/61436-regularizend
% [2] Gridfit. Matlab File Exchange open source. 
% https://www.mathworks.com/matlabcentral/fileexchange/8998-surface-fitting-using-gridfit
% ----------------------------------------------
% Author: Jin Yang.
% Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
% Last time updated: 2020.12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialization
h = DICmesh.elementMinSize;
coordinatesFEM = DICmesh.coordinatesFEM;
FilterSizeInput = DICpara.StrainFilterSize;
FilterStd = DICpara.StrainFilterStd; 
F = full(F); 
try smoothness = DICpara.StrainSmoothness; 
catch smoothness = 1e-5;
end


%% prompt = 'Do you want to smooth displacement? (0-yes; 1-no)';
DoYouWantToSmoothOnceMore = 0; % DoYouWantToSmoothOnceMore = input(prompt);
if DoYouWantToSmoothOnceMore == 0  
    if isempty(FilterStd) == 1
        prompt = 'Choose filter standard deviation(0-default): ';
        FilterStd = input(prompt);
        if FilterStd == 0
            FilterStd = 0.5; 
        end
    else
        if FilterStd == 0
            FilterStd = 0.5;
        end
    end
    if isempty(FilterSizeInput) == 1
        prompt = 'Choose Gaussian filter size(0-default): ';
        FilterSizeInput = input(prompt);
        if FilterSizeInput == 0
            FilterSizeInput = 2*ceil(2*FilterStd)+1; 
        end
    else
        if FilterSizeInput == 0
            FilterSizeInput = 2*ceil(2*FilterStd)+1;
        end
    end
end

SmoothTimes = 1;
while (DoYouWantToSmoothOnceMore==0)
    Coordxnodes = [min(coordinatesFEM(:,1)):h:max(coordinatesFEM(:,1))]'; 
    Coordynodes = [min(coordinatesFEM(:,2)):h:max(coordinatesFEM(:,2))]';
    % Iblur_11 = gridfit(coordinatesFEM(:,1), coordinatesFEM(:,2), F(1:4:end),Coordxnodes,Coordynodes,'regularizer','springs'); 
    % Iblur_11=Iblur_11';
    % Iblur_22 = gridfit(coordinatesFEM(:,1), coordinatesFEM(:,2), F(4:4:end),Coordxnodes,Coordynodes,'regularizer','springs');  
    % Iblur_22=Iblur_22';
    % Iblur_21 = gridfit(coordinatesFEM(:,1), coordinatesFEM(:,2), F(2:4:end),Coordxnodes,Coordynodes,'regularizer','springs'); 
    % Iblur_21=Iblur_21';
    % Iblur_12 = gridfit(coordinatesFEM(:,1), coordinatesFEM(:,2), F(3:4:end),Coordxnodes,Coordynodes,'regularizer','springs');
    % Iblur_12=Iblur_12';
    
    Iblur_11 = regularizeNd([coordinatesFEM(:,1), coordinatesFEM(:,2)],F(1:4:end),{Coordxnodes,Coordynodes},smoothness);
    Iblur_22 = regularizeNd([coordinatesFEM(:,1), coordinatesFEM(:,2)],F(4:4:end),{Coordxnodes,Coordynodes},smoothness);
    Iblur_21 = regularizeNd([coordinatesFEM(:,1), coordinatesFEM(:,2)],F(2:4:end),{Coordxnodes,Coordynodes},smoothness);
    Iblur_12 = regularizeNd([coordinatesFEM(:,1), coordinatesFEM(:,2)],F(3:4:end),{Coordxnodes,Coordynodes},smoothness);
    
    % -------------------------------------------------------
    imageFilter=fspecial('gaussian',FilterSizeInput,FilterStd);
    Iblur_1 = nanconv(Iblur_11,imageFilter,'edge','nanout');
    Iblur_4 = nanconv(Iblur_22,imageFilter,'edge','nanout');
    Iblur_2 = nanconv(Iblur_21,imageFilter,'edge','nanout');
    Iblur_3 = nanconv(Iblur_12,imageFilter,'edge','nanout');
    
    for tempi = 1:size(coordinatesFEM,1)
        [row1,~] = find(Coordxnodes==coordinatesFEM(tempi,1));
        [row2,~] = find(Coordynodes==coordinatesFEM(tempi,2));
        F(4*tempi-3) = Iblur_1(row1,row2);
        F(4*tempi)   = Iblur_4(row1,row2);
        F(4*tempi-2) = Iblur_2(row1,row2);
        F(4*tempi-1) = Iblur_3(row1,row2);
    end
     
    % prompt = 'Do you want to smooth displacement once more? (0-yes; 1-no)';
    % DoYouWantToSmoothOnceMore = input(prompt); 
    SmoothTimes = SmoothTimes+1;
    if SmoothTimes > 2
        DoYouWantToSmoothOnceMore = 1;
    end
    
end

