function [StrainNodalPt,StrainGaussPt,CoordGaussPt] = funGlobalNodalStrainQuadtree(...
                                           DICmesh,U,GaussPtOrder,waitBarDisplayOrNot)
%FUNGLOBALNODALSTRAINQUADTREE:  to compute strain fields by the FE-method
%   [StrainNodalPt,StrainGaussPt,CoordGaussPt] = funGlobalNodalStrainQuadtree( ...
%                                             DICmesh,U,GaussPtOrder,waitBarDisplayOrNot)                 
% ----------------------------------------------
%
%	INPUT: DICmesh             DIC FE Q4 mesh: coordinatesFEM, elementsFEM
%          U                   Disp vector: U = [Ux_node1, Uy_node1, ... , Ux_nodeN, Uy_nodeN]';
%          GaussPtOrder        Gauss point order used in FE Q4 element 
%          udual,vdual         Dual variables
%          waitBarDisplayOrNot Display a waitbar or not
%           
%   OUTPUT: StrainNodalPt      Solved strains at the nodal points
%           StrainGaussPt      Solved strains at the Gauss points
%           CoordGaussPt       Coordinates of Gauss points
%
% ----------------------------------------------
% Reference
% [1] RegularizeNd. Matlab File Exchange open source. 
% https://www.mathworks.com/matlabcentral/fileexchange/61436-regularizend
% [2] Gridfit. Matlab File Exchange open source. 
% https://www.mathworks.com/matlabcentral/fileexchange/8998-surface-fitting-using-gridfit
% [3] Rbfinterp. Matlab File Exchange open source.
% https://www.mathworks.com/matlabcentral/fileexchange/10056-scattered-data-interpolation-and-approximation-using-radial-base-functions
% ----------------------------------------------
% Author: Jin Yang.  
% Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
% Last time updated: 2018.03, 2020.12 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialization
coordinatesFEM = DICmesh.coordinatesFEM; 
elementsFEM = DICmesh.elementsFEM; elementsFEM(:,5:8) = 0;
h = DICmesh.elementMinSize;

DIM = 2; NodesNumPerEle = 4;
SizeCoords = size(coordinatesFEM,1);
smoothness = 0; % an empirical value; if apply some smoothness, use "1e-3" by default


%% ====== Initialize strain Gauss points ======
U = [U;zeros(DIM*NodesNumPerEle,1)]; % For zero entries in elements
StrainGaussPt = zeros(4*size(elementsFEM,1),4); CoordGaussPt = zeros(4*size(elementsFEM,1),2);

% ====== Gaussian quadrature parameter ======
switch GaussPtOrder
    case 2 % 2*2 Gauss points 
        gqpt1 = -1/sqrt(3); gqpt2 = 1/sqrt(3); gqpt = [gqpt1,gqpt2]; 
        gqwt1 = 1; gqwt2 = 1; gqwt = [gqwt1,gqwt2];
    case 3 % 3*3 Gauss points 
        gqpt1 = 0; gqpt2 = sqrt(3/5); gqpt3 = -sqrt(3/5); gqpt = [gqpt1,gqpt2,gqpt3];
        gqwt1 = 8/9; gqwt2 = 5/9; gqwt3 = 5/9; gqwt = [gqwt1,gqwt2,gqwt3];
    case 4 % 4*4 Gauss points 
        gqpt1 = 0.339981; gqpt2 = -0.339981; gqpt3 = 0.861136; gqpt4 = -0.861136;  
        gqwt1 = 0.652145; gqwt2 = 0.652145; gqwt3 = 0.347855; gqwt4 = 0.347855;
        gqpt = [gqpt1,gqpt2,gqpt3,gqpt4]; gqwt = [gqwt1,gqwt2,gqwt3,gqwt4];
    case 5 % 5*5 Gauss points 
        gqpt1 = 0; gqpt2 = 0.538469; gqpt3 = -0.538469; gqpt4 = 0.90618; gqpt5 = -0.90618;
        gqwt1 = 0.568889; gqwt2 = 0.478629; gqwt3 = 0.478629; gqwt4 = 0.236927; gqwt5 = 0.236927;
        gqpt = [gqpt1,gqpt2,gqpt3,gqpt4,gqpt5]; gqwt = [gqwt1,gqwt2,gqwt3,gqwt4,gqwt5];
    otherwise
        disp('Not implemented yet!')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
if waitBarDisplayOrNot == 0
   hbar = waitbar(0, ['Update deformation gradient tensor F after solving Subproblem 2.']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for eleInd = 1:size(elementsFEM,1) % eleInd is the element index
    
    if waitBarDisplayOrNot == 0
        waitbar(eleInd/size(elementsFEM,1));
    end
    
    StrainWithinEachElementGausspoint = zeros(GaussPtOrder^2,4);
    
    % ----- Find four corner points -------
    pt1x = coordinatesFEM(elementsFEM(eleInd,1),1); pt1y = coordinatesFEM(elementsFEM(eleInd,1),2);
    pt2x = coordinatesFEM(elementsFEM(eleInd,2),1); pt2y = coordinatesFEM(elementsFEM(eleInd,2),2);
    pt3x = coordinatesFEM(elementsFEM(eleInd,3),1); pt3y = coordinatesFEM(elementsFEM(eleInd,3),2);
    pt4x = coordinatesFEM(elementsFEM(eleInd,4),1); pt4y = coordinatesFEM(elementsFEM(eleInd,4),2);
    
    % ------ Find mid points 5/6/7/8 -------
    if elementsFEM(eleInd,5) ~= 0
        pt5x = coordinatesFEM(elementsFEM(eleInd,5),1); pt5y = coordinatesFEM(elementsFEM(eleInd,5),2);
    else, pt5x = 0; pt5y = 0; 
    end
    if elementsFEM(eleInd,6) ~= 0
        pt6x = coordinatesFEM(elementsFEM(eleInd,6),1); pt6y = coordinatesFEM(elementsFEM(eleInd,6),2);
    else, pt6x = 0; pt6y = 0; 
    end
    if elementsFEM(eleInd,7) ~= 0
        pt7x = coordinatesFEM(elementsFEM(eleInd,7),1); pt7y = coordinatesFEM(elementsFEM(eleInd,7),2);
    else, pt7x = 0; pt7y = 0; 
    end
    if elementsFEM(eleInd,8) ~= 0
        pt8x = coordinatesFEM(elementsFEM(eleInd,8),1); pt8y = coordinatesFEM(elementsFEM(eleInd,8),2);
    else, pt8x = 0; pt8y = 0; 
    end

    % ------ Calculate ksi and eta --------
    lMatrix = [ pt1x*pt1y pt1x pt1y 1;
                pt2x*pt2y pt2x pt2y 1;
                pt3x*pt3y pt3x pt3y 1;
                pt4x*pt4y pt4x pt4y 1 ];

    % ------ Find the element nodal indices ------
    lb = [-1;1;1;-1]; l = linsolve(lMatrix,lb);
    mb = [-1;-1;1;1]; m = linsolve(lMatrix,mb);

    % ------ Find the element nodal indices ------
    tempIndexU = 2*elementsFEM(eleInd,[1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8]);
    tempIndexU(1:2:end) = tempIndexU(1:2:end)-1;   
    
    % We don't want temp <= 0, instead, put them to the end
    for tempi = 5:8
        if elementsFEM(eleInd,tempi) == 0
            tempIndexU(2*tempi-1) = 2*(SizeCoords+1)-1;    
            tempIndexU(2*tempi)   = 2*(SizeCoords+1);
        end
    end
     
    % ------ Set Gauss points ------
    pt1 = elementsFEM(eleInd,1); pt2 = elementsFEM(eleInd,3);
    ptOfx = zeros(length(gqwt),1); ptOfy = zeros(length(gqwt),1);
    for tempi = 1:length(gqwt)
       ptOfx(tempi) = gqpt(tempi)*0.5*(coordinatesFEM(pt2,1)-coordinatesFEM(pt1,1))+0.5*(coordinatesFEM(pt2,1)+coordinatesFEM(pt1,1));
       ptOfy(tempi) = gqpt(tempi)*0.5*(coordinatesFEM(pt2,2)-coordinatesFEM(pt1,2))+0.5*(coordinatesFEM(pt2,2)+coordinatesFEM(pt1,2));
    end
    
    % ------ Start computing strains at Gauss points ------
    for tempi = 1:length(ptOfx)
       for tempj = 1:length(ptOfy)
     
            % ------ Calculate ksi and eta ------
            ksi = l(1)*ptOfx(tempi)*ptOfy(tempj) + l(2)*ptOfx(tempi) + l(3)*ptOfy(tempj) + l(4);
            eta = m(1)*ptOfx(tempi)*ptOfy(tempj) + m(2)*ptOfx(tempi) + m(3)*ptOfy(tempj) + m(4);
             
            % ------ Calculate N ------
            deltaPt5 = 0; deltaPt6 = 0; deltaPt7 = 0; deltaPt8 = 0;
            if elementsFEM(eleInd,5) ~= 0; deltaPt5 = 1; end
            if elementsFEM(eleInd,6) ~= 0; deltaPt6 = 1; end
            if elementsFEM(eleInd,7) ~= 0; deltaPt7 = 1; end
            if elementsFEM(eleInd,8) ~= 0; deltaPt8 = 1; end

            N5 = deltaPt5*0.5*(1+ksi)*(1-abs(eta));
            N6 = deltaPt6*0.5*(1+eta)*(1-abs(ksi));
            N7 = deltaPt7*0.5*(1-ksi)*(1-abs(eta));
            N8 = deltaPt8*0.5*(1-eta)*(1-abs(ksi));

            N1 = (1-ksi)*(1-eta)*0.25 - 0.5*(N7+N8);
            N2 = (1+ksi)*(1-eta)*0.25 - 0.5*(N8+N5);
            N3 = (1+ksi)*(1+eta)*0.25 - 0.5*(N5+N6);
            N4 = (1-ksi)*(1+eta)*0.25 - 0.5*(N6+N7);

            % ------ Build J matrix ------
            % Jacobian matrix J doesn't change with added enriched functions  
            J11 = funDN1Dksi(ksi,eta,deltaPt7,deltaPt8)*pt1x + funDN2Dksi(ksi,eta,deltaPt8,deltaPt5)*pt2x + ...
                  funDN3Dksi(ksi,eta,deltaPt5,deltaPt6)*pt3x + funDN4Dksi(ksi,eta,deltaPt6,deltaPt7)*pt4x + ...
                  funDN5Dksi(ksi,eta,deltaPt5)*pt5x + funDN6Dksi(ksi,eta,deltaPt6)*pt6x + ...
                  funDN7Dksi(ksi,eta,deltaPt7)*pt7x + funDN8Dksi(ksi,eta,deltaPt8)*pt8x;
            J12 = funDN1Dksi(ksi,eta,deltaPt7,deltaPt8)*pt1y + funDN2Dksi(ksi,eta,deltaPt8,deltaPt5)*pt2y + ...
                  funDN3Dksi(ksi,eta,deltaPt5,deltaPt6)*pt3y + funDN4Dksi(ksi,eta,deltaPt6,deltaPt7)*pt4y + ...
                  funDN5Dksi(ksi,eta,deltaPt5)*pt5y + funDN6Dksi(ksi,eta,deltaPt6)*pt6y + ...
                  funDN7Dksi(ksi,eta,deltaPt7)*pt7y + funDN8Dksi(ksi,eta,deltaPt8)*pt8y;
            J21 = funDN1Deta(ksi,eta,deltaPt7,deltaPt8)*pt1x + funDN2Deta(ksi,eta,deltaPt8,deltaPt5)*pt2x + ...
                  funDN3Deta(ksi,eta,deltaPt5,deltaPt6)*pt3x + funDN4Deta(ksi,eta,deltaPt6,deltaPt7)*pt4x + ...
                  funDN5Deta(ksi,eta,deltaPt5)*pt5x + funDN6Deta(ksi,eta,deltaPt6)*pt6x + ...
                  funDN7Deta(ksi,eta,deltaPt7)*pt7x + funDN8Deta(ksi,eta,deltaPt8)*pt8x;
            J22 = funDN1Deta(ksi,eta,deltaPt7,deltaPt8)*pt1y + funDN2Deta(ksi,eta,deltaPt8,deltaPt5)*pt2y + ...
                  funDN3Deta(ksi,eta,deltaPt5,deltaPt6)*pt3y + funDN4Deta(ksi,eta,deltaPt6,deltaPt7)*pt4y + ...
                  funDN5Deta(ksi,eta,deltaPt5)*pt5y + funDN6Deta(ksi,eta,deltaPt6)*pt6y + ...
                  funDN7Deta(ksi,eta,deltaPt7)*pt7y + funDN8Deta(ksi,eta,deltaPt8)*pt8y;

            J = [J11 J12; J21 J22];
            Jacobian = det(J);
            InvJ = 1/Jacobian*[J22 -J12; -J21 J11];
            
            % ------ Compute DN matrix ------
            DN = [InvJ zeros(2,2); zeros(2,2) InvJ] * ...
                    [funDN1Dksi(ksi,eta,deltaPt7,deltaPt8) 0 funDN2Dksi(ksi,eta,deltaPt8,deltaPt5) 0 ...
                    funDN3Dksi(ksi,eta,deltaPt5,deltaPt6) 0 funDN4Dksi(ksi,eta,deltaPt6,deltaPt7) 0 ...
                    funDN5Dksi(ksi,eta,deltaPt5) 0 funDN6Dksi(ksi,eta,deltaPt6) 0 ...
                    funDN7Dksi(ksi,eta,deltaPt7) 0 funDN8Dksi(ksi,eta,deltaPt8) 0;
                    funDN1Deta(ksi,eta,deltaPt7,deltaPt8) 0 funDN2Deta(ksi,eta,deltaPt8,deltaPt5) 0 ...
                    funDN3Deta(ksi,eta,deltaPt5,deltaPt6) 0 funDN4Deta(ksi,eta,deltaPt6,deltaPt7) 0 ...
                    funDN5Deta(ksi,eta,deltaPt5) 0 funDN6Deta(ksi,eta,deltaPt6) 0 ...
                    funDN7Deta(ksi,eta,deltaPt7) 0 funDN8Deta(ksi,eta,deltaPt8) 0;
                    0 funDN1Dksi(ksi,eta,deltaPt7,deltaPt8) 0 funDN2Dksi(ksi,eta,deltaPt8,deltaPt5) ...
                    0 funDN3Dksi(ksi,eta,deltaPt5,deltaPt6) 0 funDN4Dksi(ksi,eta,deltaPt6,deltaPt7) ...
                    0 funDN5Dksi(ksi,eta,deltaPt5) 0 funDN6Dksi(ksi,eta,deltaPt6) ...
                    0 funDN7Dksi(ksi,eta,deltaPt7) 0 funDN8Dksi(ksi,eta,deltaPt8);
                    0 funDN1Deta(ksi,eta,deltaPt7,deltaPt8) 0 funDN2Deta(ksi,eta,deltaPt8,deltaPt5) ...
                    0 funDN3Deta(ksi,eta,deltaPt5,deltaPt6) 0 funDN4Deta(ksi,eta,deltaPt6,deltaPt7) ...
                    0 funDN5Deta(ksi,eta,deltaPt5) 0 funDN6Deta(ksi,eta,deltaPt6) ...
                    0 funDN7Deta(ksi,eta,deltaPt7) 0 funDN8Deta(ksi,eta,deltaPt8);];
             
            StrainWithinEachElementGausspoint(length(gqpt)*(tempi-1)+tempj,1:4) = DN*U(tempIndexU);
            
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    switch GaussPtOrder 
        case 2 % 2*2 Gauss point 
            StrainGaussPt((length(gqpt)^2*(eleInd-1)+1):length(gqpt)^2*eleInd,1:4) = StrainWithinEachElementGausspoint;
            CoordGaussPt((length(gqpt)^2*(eleInd-1)+1):length(gqpt)^2*eleInd,1:2) = ...
                    [ptOfx(1),ptOfy(1);ptOfx(1),ptOfy(2); 
                     ptOfx(2),ptOfy(1);ptOfx(2),ptOfy(2)];
        case 3 % 3*3 Gauss point 
            StrainGaussPt((length(gqpt)^2*(eleInd-1)+1):length(gqpt)^2*eleInd,1:4) = StrainWithinEachElementGausspoint;
            CoordGaussPt((length(gqpt)^2*(eleInd-1)+1):length(gqpt)^2*eleInd,1:2) = ...
                    [ptOfx(1),ptOfy(1);ptOfx(1),ptOfy(2);ptOfx(1),ptOfy(3);
                     ptOfx(2),ptOfy(1);ptOfx(2),ptOfy(2);ptOfx(2),ptOfy(3);
                     ptOfx(3),ptOfy(1);ptOfx(3),ptOfy(2);ptOfx(3),ptOfy(3)];
        case 4 % 4*4 Gauss point
            StrainGaussPt((length(gqpt)^2*(eleInd-1)+1):length(gqpt)^2*eleInd,1:4) = StrainWithinEachElementGausspoint;
            CoordGaussPt((length(gqpt)^2*(eleInd-1)+1):length(gqpt)^2*eleInd,1:2) = ...
                    [ptOfx(1),ptOfy(1);ptOfx(1),ptOfy(2);ptOfx(1),ptOfy(3);ptOfx(1),ptOfy(4);
                     ptOfx(2),ptOfy(1);ptOfx(2),ptOfy(2);ptOfx(2),ptOfy(3);ptOfx(2),ptOfy(4);
                     ptOfx(3),ptOfy(1);ptOfx(3),ptOfy(2);ptOfx(3),ptOfy(3);ptOfx(3),ptOfy(4);
                     ptOfx(4),ptOfy(1);ptOfx(4),ptOfy(2);ptOfx(4),ptOfy(3);ptOfx(4),ptOfy(4)];
        case 5 % 5*5 Gauss point
            StrainGaussPt((length(gqpt)^2*(eleInd-1)+1):length(gqpt)^2*eleInd,1:4) = StrainWithinEachElementGausspoint;
            CoordGaussPt((length(gqpt)^2*(eleInd-1)+1):length(gqpt)^2*eleInd,1:2) = ...
               [ptOfx(1),ptOfy(1);ptOfx(1),ptOfy(2);ptOfx(1),ptOfy(3);ptOfx(1),ptOfy(4);ptOfx(1),ptOfy(5);
                ptOfx(2),ptOfy(1);ptOfx(2),ptOfy(2);ptOfx(2),ptOfy(3);ptOfx(2),ptOfy(4);ptOfx(2),ptOfy(5);
                ptOfx(3),ptOfy(1);ptOfx(3),ptOfy(2);ptOfx(3),ptOfy(3);ptOfx(3),ptOfy(4);ptOfx(3),ptOfy(5);
                ptOfx(4),ptOfy(1);ptOfx(4),ptOfy(2);ptOfx(4),ptOfy(3);ptOfx(4),ptOfy(4);ptOfx(4),ptOfy(5);
                ptOfx(5),ptOfy(1);ptOfx(5),ptOfy(2);ptOfx(5),ptOfy(3);ptOfx(5),ptOfy(4);ptOfx(5),ptOfy(5);];
        otherwise
            disp('Not implemented yet!')
    end
    %%%%%%%%%%%%%%%%%%%%%%%% End of Comment %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

if waitBarDisplayOrNot==0, close(hbar); end
 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if smoothness == 0 % No smoothness
    
    %     StrainNodalPt = zeros(4*size(coordinatesFEM,1),1);
    %     
    %     dilatedI = ( imgaussfilt(double(DICpara.ImgRefMask),1) );
    %     dilatedI = logical( dilatedI > 0.01);
    %     cc = bwconncomp(dilatedI,8);
    %     indPxAll = sub2ind( DICpara.ImgSize, round(coordinatesFEM(:,1)), round(coordinatesFEM(:,2)) );
    %     indPxNotNanAll = sub2ind( DICpara.ImgSize, round(CoordGaussPt(:,1)), round(CoordGaussPt(:,2)) );
    %     stats = regionprops(cc,'Area','PixelList');
    %     for tempi = 1:length(stats)
    %     
    %         if stats(tempi).Area > 20
    %             
    %             %%%%% Find those nodes %%%%%
    %             indPxtempi = sub2ind( DICpara.ImgSize, stats(tempi).PixelList(:,2), stats(tempi).PixelList(:,1) );
    %             Lia = ismember(indPxAll,indPxtempi); [LiaList,~] = find(Lia==1);
    %             Lib = ismember(indPxNotNanAll,indPxtempi); [LibList,~] = find(Lib==1);
    %             
    %             % ------ F11 ------
    %             op1 = rbfcreate( [CoordGaussPt( (LibList),1:2)]',[StrainGaussPt(LibList,1)]','RBFFunction', 'thinplate'); rbfcheck(op1);
    %             fi1 = rbfinterp( [coordinatesFEM(LiaList,1:2)]', op1);
    %             StrainNodalPt(4*LiaList-3) = fi1(:);
    %             % ------ F21 ------
    %             op1 = rbfcreate( [CoordGaussPt( (LibList),1:2)]',[StrainGaussPt(LibList,2)]','RBFFunction', 'thinplate'); rbfcheck(op1);
    %             fi1 = rbfinterp( [coordinatesFEM(LiaList,1:2)]', op1);
    %             StrainNodalPt(4*LiaList-2) = fi1(:);
    %             % ------ F12 ------
    %             op1 = rbfcreate( [CoordGaussPt( (LibList),1:2)]',[StrainGaussPt(LibList,3)]','RBFFunction', 'thinplate'); rbfcheck(op1);
    %             fi1 = rbfinterp( [coordinatesFEM(LiaList,1:2)]', op1);
    %             StrainNodalPt(4*LiaList-1) = fi1(:);
    %             % ------ F22 ------
    %             op1 = rbfcreate( [CoordGaussPt( (LibList),1:2)]',[StrainGaussPt(LibList,4)]','RBFFunction', 'thinplate'); rbfcheck(op1);
    %             fi1 = rbfinterp( [coordinatesFEM(LiaList,1:2)]', op1);
    %             StrainNodalPt(4*LiaList) = fi1(:);
    %             
    %         end
    %     end
    
    Ftemp = scatteredInterpolant(CoordGaussPt(:,1), CoordGaussPt(:,2), StrainGaussPt(:,1),'natural','nearest');
    F11 = Ftemp(coordinatesFEM(:,1),coordinatesFEM(:,2));
    Ftemp = scatteredInterpolant(CoordGaussPt(:,1), CoordGaussPt(:,2), StrainGaussPt(:,2),'natural','nearest');
    F21 = Ftemp(coordinatesFEM(:,1),coordinatesFEM(:,2));
    Ftemp = scatteredInterpolant(CoordGaussPt(:,1), CoordGaussPt(:,2), StrainGaussPt(:,3),'natural','nearest');
    F12 = Ftemp(coordinatesFEM(:,1),coordinatesFEM(:,2));
    Ftemp = scatteredInterpolant(CoordGaussPt(:,1), CoordGaussPt(:,2), StrainGaussPt(:,4),'natural','nearest');
    F22 = Ftemp(coordinatesFEM(:,1),coordinatesFEM(:,2));
    
    StrainNodalPt = [F11(:),F21(:),F12(:),F22(:)]'; StrainNodalPt = StrainNodalPt(:);
    
else %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % ====== Use gridfit and Gauss points to fit strain field ======
    Coordxnodes = [min(coordinatesFEM(:,1)):h:max(coordinatesFEM(:,1))]';
    Coordynodes = [min(coordinatesFEM(:,2)):h:max(coordinatesFEM(:,2))]';
    disp('** This may take a while when quadtree minimum element size is small. Please wait. ')
    
    Iblur_10 = regularizeNd([CoordGaussPt(:,1),CoordGaussPt(:,2)],StrainGaussPt(:,1),{Coordxnodes,Coordynodes},smoothness);
    disp('|------> 25%                   | Please keep waiting! ')
    Iblur_20 = regularizeNd([CoordGaussPt(:,1),CoordGaussPt(:,2)],StrainGaussPt(:,2),{Coordxnodes,Coordynodes},smoothness);
    disp('|-------------> 50%            | Half way! Do not give up! ')
    Iblur_30 = regularizeNd([CoordGaussPt(:,1),CoordGaussPt(:,2)],StrainGaussPt(:,3),{Coordxnodes,Coordynodes},smoothness);
    disp('|--------------------> 75%     | Almost there! ')
    Iblur_40 = regularizeNd([CoordGaussPt(:,1),CoordGaussPt(:,2)],StrainGaussPt(:,4),{Coordxnodes,Coordynodes},smoothness);
    disp('|----------------------------->| Hooray!!! ')
     
    % Iblur_10 = gridfit(CoordGaussPt(:,1),CoordGaussPt(:,2),StrainGaussPt(:,1),Coordxnodes,Coordynodes ); Iblur_10=Iblur_10';
    % disp('|------> 25%                   | Please keep waiting! ')
    % Iblur_20 = gridfit(CoordGaussPt(:,1),CoordGaussPt(:,2),StrainGaussPt(:,2),Coordxnodes,Coordynodes ); Iblur_20=Iblur_20';
    % disp('|-------------> 50%            | Half way! Do not give up! ')
    % Iblur_30 = gridfit(CoordGaussPt(:,1),CoordGaussPt(:,2),StrainGaussPt(:,3),Coordxnodes,Coordynodes ); Iblur_30=Iblur_30';
    % disp('|--------------------> 75%     | Almost there! ')
    % Iblur_40 = gridfit(CoordGaussPt(:,1),CoordGaussPt(:,2),StrainGaussPt(:,4),Coordxnodes,Coordynodes ); Iblur_40=Iblur_40';
    % disp('|----------------------------->| Hoooooray!!! ')
     
    StrainNodalPt = 0*[U(1:end-DIM*NodesNumPerEle); U(1:end-DIM*NodesNumPerEle)];
    for tempi = 1:size(coordinatesFEM,1)
        [row1,~] = find(Coordxnodes == coordinatesFEM(tempi,1));
        [row2,~] = find(Coordynodes == coordinatesFEM(tempi,2));
        StrainNodalPt(4*tempi-3) = Iblur_10(row1,row2);
        StrainNodalPt(4*tempi-2) = Iblur_20(row1,row2);
        StrainNodalPt(4*tempi-1) = Iblur_30(row1,row2);
        StrainNodalPt(4*tempi-0) = Iblur_40(row1,row2);
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end



%% ========= subroutines for  FEM Q4 shape function derivatives ========
function DN5Dksi=funDN5Dksi(ksi,eta,deltaPt5)
DN5Dksi = deltaPt5*0.5*(1-abs(eta)) ;
end
function DN5Deta=funDN5Deta(ksi,eta,deltaPt5)
DN5Deta = deltaPt5*0.5*(1+ksi)*sign(-eta);
end
function DN6Dksi=funDN6Dksi(ksi,eta,deltaPt6)
DN6Dksi = deltaPt6*0.5*(1+eta)*sign(-ksi);
end
function DN6Deta=funDN6Deta(ksi,eta,deltaPt6)
DN6Deta = deltaPt6*0.5*(1-abs(ksi));
end
function DN7Dksi=funDN7Dksi(ksi,eta,deltaPt7)
DN7Dksi = deltaPt7*0.5*(-1)*(1-abs(eta));
end
function DN7Deta=funDN7Deta(ksi,eta,deltaPt7)
DN7Deta = deltaPt7*0.5*(1-ksi)*sign(-eta);
end
function DN8Dksi=funDN8Dksi(ksi,eta,deltaPt8)
DN8Dksi = deltaPt8*0.5*(1-eta)*sign(-ksi);
end
function DN8Deta=funDN8Deta(ksi,eta,deltaPt8)
DN8Deta = deltaPt8*0.5*(-1)*(1-abs(ksi));
end
function DN1Dksi = funDN1Dksi(ksi,eta,deltaPt7,deltaPt8)
DN1Dksi = -0.25*(1-eta)-0.5*((deltaPt7*0.5*(-1)*(1-abs(eta)))+deltaPt8*0.5*(1-eta)*sign(-ksi));
end
function DN1Deta = funDN1Deta(ksi,eta,deltaPt7,deltaPt8)
DN1Deta = -0.25*(1-ksi)-0.5*(deltaPt7*0.5*(1-ksi)*sign(-eta)+deltaPt8*0.5*(-1)*(1-abs(ksi)));
end
function DN2Dksi = funDN2Dksi(ksi,eta,deltaPt8,deltaPt5)
DN2Dksi = 0.25*(1-eta)-0.5*(deltaPt8*0.5*(1-eta)*sign(-ksi)+deltaPt5*0.5*(1-abs(eta)));
end
function DN2Deta = funDN2Deta(ksi,eta,deltaPt8,deltaPt5)
DN2Deta = -0.25*(1+ksi)-0.5*(deltaPt8*0.5*(-1)*(1-abs(ksi))+deltaPt5*0.5*(1+ksi)*sign(-eta));
end
function DN3Dksi = funDN3Dksi(ksi,eta,deltaPt5,deltaPt6)
DN3Dksi = 0.25*(1+eta)-0.5*(deltaPt5*0.5*(1-abs(eta))+deltaPt6*0.5*(1+eta)*sign(-ksi));
end
function DN3Deta = funDN3Deta(ksi,eta,deltaPt5,deltaPt6)
DN3Deta = 0.25*(1+ksi)-0.5*(deltaPt5*0.5*(1+ksi)*sign(-eta)+deltaPt6*0.5*(1-abs(ksi)));
end
function DN4Dksi = funDN4Dksi(ksi,eta,deltaPt6,deltaPt7)
DN4Dksi = -0.25*(1+eta)-0.5*(deltaPt6*0.5*(1+eta)*sign(-ksi)+deltaPt7*0.5*(-1)*(1-abs(eta)));
end
function DN4Deta = funDN4Deta(ksi,eta,deltaPt6,deltaPt7)
DN4Deta = 0.25*(1-ksi)-0.5*(deltaPt6*0.5*(1-abs(ksi))+deltaPt7*0.5*(1-ksi)*sign(-eta));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 