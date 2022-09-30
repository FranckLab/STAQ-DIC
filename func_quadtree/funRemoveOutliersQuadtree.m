function [U,F] = funRemoveOutliersQuadtree(DICmesh,DICpara,U,F,varargin)
% =========================================================================
% removes outliers using the universal
% outlier test based on
%
% J. Westerweel and F. Scarano. Universal outlier detection for PIV data.
% Exp. Fluids, 39(6):1096{1100, August 2005. doi: 10.1007/s00348-005-0016-6
% -------------------------------------------------------------------------
% NOTES
% -------------------------------------------------------------------------
% needs medFilt3 and John D'Errico's inpaint_nans3
% (http://www.mathworks.com/matlabcentral/fileexchange/4551-inpaint-nans)function.
% =========================================================================

winstepsize = DICpara.winstepsize;
coordinatesFEM = DICmesh.coordinatesFEM;
elementsFEM = DICmesh.elementsFEM;

switch nargin
    case 7
        BadptRow = varargin{1}; BadptCol = varargin{2};
    otherwise
        BadptRow = []; BadptCol = [];
end

%%
close all;
Ux = U(1:2:end);
Sqx = zeros(4,size(elementsFEM,1)); Sqy = zeros(4,size(elementsFEM,1)); Sqc = zeros(4,size(elementsFEM,1));
for j = 1:size(elementsFEM,1)
    Sqx(1:4,j) = coordinatesFEM(elementsFEM(j,1:4),1);
    Sqy(1:4,j) = coordinatesFEM(elementsFEM(j,1:4),2);
    Sqc(1:4,j) = Ux(elementsFEM(j,1:4));
end
if size(elementsFEM,1)>2e4
    patch(Sqx,Sqy,Sqc,'facecolor','interp','edgecolor','none');
else
    patch(Sqx,Sqy,Sqc,'facecolor','interp');
end
view(2); axis tight; axis equal; colorbar; xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
set(gcf,'color','w'); colormap jet;
title('Click bad dispx-u points and press -Enter- after done','fontweight','normal')


%% %%%%% Remove bad x-disp bad points %%%%%
fprintf('Do you clear bad points by directly pointing x-disp bad points? (0-yes; 1-no)  \n')
prompt = 'Input here: ';
ClearBadInitialPointsOrNot = input(prompt);

while ClearBadInitialPointsOrNot == 0
    
    [row1, col1] = ginput; row = floor(col1); col = floor(row1);
    
    % Make sure it's within the image domain
    BadptCol=[row]; BadptRow=[col];
    [tempindex1] = find(BadptRow<1+max(coordinatesFEM(:,1)));  [tempindex2] = find(BadptCol<1+max(coordinatesFEM(:,2)));
    [tempindex3] = find(BadptRow>min(coordinatesFEM(:,1)));  [tempindex4] = find(BadptCol>min(coordinatesFEM(:,2)));
    tempindex1 = reshape(tempindex1,length(tempindex1),1);
    tempindex2 = reshape(tempindex2,length(tempindex2),1);
    tempindex3 = reshape(tempindex3,length(tempindex3),1);
    tempindex4 = reshape(tempindex4,length(tempindex4),1);
    tempindex = unique(intersect(tempindex4,intersect(tempindex3,intersect(tempindex1,tempindex2))));
    
    BadptRow = BadptRow(tempindex); BadptCol = BadptCol(tempindex);
    row=BadptRow; col=BadptCol;
    
    % --- Find elements include clicked bad points ---
    elementFEMCenterCoordx = 0.25* (coordinatesFEM( elementsFEM(:,1)  , 1) + coordinatesFEM( elementsFEM(:,2)  , 1) + ...
        coordinatesFEM( elementsFEM(:,3)  , 1) + coordinatesFEM( elementsFEM(:,4)  , 1));
    elementFEMCenterCoordy = 0.25* (coordinatesFEM( elementsFEM(:,1)  , 2) + coordinatesFEM( elementsFEM(:,2)  , 2) + ...
        coordinatesFEM( elementsFEM(:,3)  , 2) + coordinatesFEM( elementsFEM(:,4)  , 2));
    
    BadptCoord = [];
    for tempk = 1:length(row) % iterate for each clicking points
        
        DistElementAll = sqrt( ( elementFEMCenterCoordx-row(tempk) ).^2 + (elementFEMCenterCoordy-col(tempk)).^2 );
        DistCoordAll = sqrt( ( coordinatesFEM(:,1)-row(tempk) ).^2 + (coordinatesFEM(:,2)-col(tempk)).^2 );
        [row1,~] = find(DistElementAll < 1.42*mean(winstepsize));
        [row2,~] = find(DistCoordAll < 1.42*mean(winstepsize));
        temp1 = elementsFEM(row1,:); temp2 = row2(:);
        
        BadptCoord = [BadptCoord; unique([temp1(:);temp2(:)])]; BadptCoord=BadptCoord(:);
        
    end
    
    % Find unique BadptCoord
    BadptCoord = unique(BadptCoord); BadptCoord = setdiff(BadptCoord,[0]);
    
    % Set bad points value to be NaNs
    U(2*BadptCoord-1) = NaN; U(2*BadptCoord) = NaN;
    F(4*BadptCoord-3) = NaN; F(4*BadptCoord-2) = NaN; F(4*BadptCoord-1) = NaN; F(4*BadptCoord) = NaN;
    
    %%%%%% Fill nans %%%%%%
    nanindex = find(isnan(U(1:2:end))==1); notnanindex = setdiff([1:1:size(coordinatesFEM,1)],nanindex);
    nanindexF = find(isnan(F(1:4:end))==1); notnanindexF = setdiff([1:1:size(coordinatesFEM,1)],nanindexF);
    
    % dilatedI = ( imgaussfilt(double(Df.ImgRefMask),0.5) );
    % dilatedI = logical( dilatedI > 0.01); % figure, imshow(dilatedI)
    dilatedI = logical(DICpara.ImgRefMask);
    cc = bwconncomp(dilatedI,8);
    indPxAll = sub2ind( DICpara.ImgSize, round(coordinatesFEM(:,1)),round(coordinatesFEM(:,2)));
    indPxNotNanAll = sub2ind( DICpara.ImgSize, round(coordinatesFEM(notnanindex,1)), round(coordinatesFEM(notnanindex,2)) );
    stats = regionprops(cc,'Area','PixelList');
    for tempi = 1:length(stats)
        
        if stats(tempi).Area > 20
            
            %%%%% Find those nodes %%%%%
            indPxtempi = sub2ind( DICpara.ImgSize, stats(tempi).PixelList(:,2), stats(tempi).PixelList(:,1) );
            Lia = ismember(indPxAll,indPxtempi); [LiaList,~] = find(Lia==1);
            Lib = ismember(indPxNotNanAll,indPxtempi); [LibList,~] = find(Lib==1);

            %%%%% RBF (Radial basis function) works better than "scatteredInterpolant" %%%%%
            % ------ Disp u ------
            op1 = rbfcreate( [coordinatesFEM(notnanindex(LibList),1:2)]',[U(2*notnanindex(LibList)-1)]','RBFFunction', 'thinplate'); rbfcheck(op1);
            fi1 = rbfinterp( [coordinatesFEM(LiaList,1:2)]', op1);
            U(2*LiaList-1) = fi1(:);
            % figure, plot3(coordinatesFEM(notnanindex(LibList),1),coordinatesFEM(notnanindex(LibList),2),U(2*notnanindex(LibList)-1),'.')
            % hold on; plot3(coordinatesFEM(LiaList,1),coordinatesFEM(LiaList,2),fi1,'.')

            % ------ Disp v ------
            op1 = rbfcreate( [coordinatesFEM(notnanindex(LibList),1:2)]',[U(2*notnanindex(LibList) )]','RBFFunction', 'thinplate'); rbfcheck(op1);
            fi1 = rbfinterp( [coordinatesFEM(LiaList,1:2)]', op1);
            U(2*LiaList ) = fi1(:);
            % ------ F11 ------
            op1 = rbfcreate( [coordinatesFEM(notnanindex(LibList),1:2)]',[F(4*notnanindex(LibList)-3)]','RBFFunction', 'thinplate'); rbfcheck(op1);
            fi1 = rbfinterp( [coordinatesFEM(LiaList,1:2)]', op1);
            F(4*LiaList-3) = fi1(:);
            % ------ F21 ------
            op1 = rbfcreate( [coordinatesFEM(notnanindex(LibList),1:2)]',[F(4*notnanindex(LibList)-2)]','RBFFunction', 'thinplate'); rbfcheck(op1);
            fi1 = rbfinterp( [coordinatesFEM(LiaList,1:2)]', op1);
            F(4*LiaList-2) = fi1(:);
            % ------ F12 ------
            op1 = rbfcreate( [coordinatesFEM(notnanindex(LibList),1:2)]',[F(4*notnanindex(LibList)-1)]','RBFFunction', 'thinplate'); rbfcheck(op1);
            fi1 = rbfinterp( [coordinatesFEM(LiaList,1:2)]', op1);
            F(4*LiaList-1) = fi1(:);
            % ------ F22 ------
            op1 = rbfcreate( [coordinatesFEM(notnanindex(LibList),1:2)]',[F(4*notnanindex(LibList)-0)]','RBFFunction', 'thinplate'); rbfcheck(op1);
            fi1 = rbfinterp( [coordinatesFEM(LiaList,1:2)]', op1);
            F(4*LiaList) = fi1(:);
        
        end
    end
    %%%%% Replace these old codes by applying RBF interpolation %%%%%
    %     % ------ inpaint nans using gridfit or scatteredInterpolant ------
    %     Coordxnodes = unique(coordinatesFEM(:,1)); Coordynodes = unique(coordinatesFEM(:,2));
    %     nanindex = find(isnan(U(1:2:end))==1); notnanindex = setdiff([1:1:size(coordinatesFEM,1)],nanindex);
    %     
    %     %[u1temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), U(2*notnanindex-1),Coordxnodes,Coordynodes,'regularizer','springs'); u1temp = u1temp';
    %     %[v1temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), U(2*notnanindex),Coordxnodes,Coordynodes,'regularizer','springs'); v1temp = v1temp';
    %     %[F11temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-3),Coordxnodes,Coordynodes,'regularizer','springs'); F11temp = F11temp';
    %     %[F21temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-2),Coordxnodes,Coordynodes,'regularizer','springs'); F21temp = F21temp';
    %     %[F12temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-1),Coordxnodes,Coordynodes,'regularizer','springs'); F12temp = F12temp';
    %     %[F22temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-0),Coordxnodes,Coordynodes,'regularizer','springs'); F22temp = F22temp';
    %     
    %     %[CoordxnodesGrid,CoordynodesGrid] = ndgrid(Coordxnodes,Coordynodes);
    %     Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), U(2*notnanindex-1),'natural','nearest');
    %     U1 = Ftemp(coordinatesFEM(:,1),coordinatesFEM(:,2));
    %     %[u1temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %u1temp = u1temp';
    %     Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), U(2*notnanindex),'natural','nearest');
    %     V1 = Ftemp(coordinatesFEM(:,1),coordinatesFEM(:,2));
    %     %[v1temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %v1temp = v1temp';
    %     Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-3),'natural','nearest');
    %     F11 = Ftemp(coordinatesFEM(:,1),coordinatesFEM(:,2));
    %     %[F11temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %F11temp = F11temp';
    %     Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-2),'natural','nearest');
    %     F21 = Ftemp(coordinatesFEM(:,1),coordinatesFEM(:,2));
    %     %[F21temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %F21temp = F21temp';
    %     Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-1),'natural','nearest');
    %     F12 = Ftemp(coordinatesFEM(:,1),coordinatesFEM(:,2));
    %     %[F12temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %F12temp = F12temp';
    %     Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-0),'natural','nearest');
    %     F22 = Ftemp(coordinatesFEM(:,1),coordinatesFEM(:,2));
    %     %[F22temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %F22temp = F22temp';
    %     
    %     %     for tempi = 1:size(coordinatesFEM,1)
    %     %         [row1,col1] = find(Coordxnodes==coordinatesFEM(tempi,1));
    %     %         [row2,col2] = find(Coordynodes==coordinatesFEM(tempi,2));
    %     %         U(2*tempi-1) = u1temp(row1,row2);
    %     %         U(2*tempi)   = v1temp(row1,row2);
    %     %         F(4*tempi-3) = F11temp(row1,row2); F(4*tempi-2) = F21temp(row1,row2);
    %     %         F(4*tempi-1) = F12temp(row1,row2); F(4*tempi) = F22temp(row1,row2);
    %     %     end
    %     U = [U1(:),V1(:)]'; U = U(:);
    %     F = [F11(:),F21(:),F12(:),F22(:)]'; F = F(:);
    
    close all;
    Ux = U(1:2:end);
    Sqx = zeros(4,size(elementsFEM,1)); Sqy = zeros(4,size(elementsFEM,1)); Sqc = zeros(4,size(elementsFEM,1));
    for j = 1:size(elementsFEM,1)
        Sqx(1:4,j) = coordinatesFEM(elementsFEM(j,1:4),1);
        Sqy(1:4,j) = coordinatesFEM(elementsFEM(j,1:4),2);
        Sqc(1:4,j) = Ux(elementsFEM(j,1:4));
    end
    if size(elementsFEM,1)>2e4
        patch(Sqx,Sqy,Sqc,'facecolor','interp','edgecolor','none');
    else
        patch(Sqx,Sqy,Sqc,'facecolor','interp');
    end
    view(2); axis tight; axis equal; colorbar; xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
    set(gcf,'color','w'); colormap jet;
    title('Click bad dispx-u points and press -Enter- after done','fontweight','normal')
    
    
    prompt = 'Do you point out more x-disp bad points? (0-yes; 1-no) Input: ';
    ClearBadInitialPointsOrNot = input(prompt);
    
end




%% %%%%% Remove bad y-disp bad points %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Have a look at integer search
% --------------------------------------
close all;
Uy = U(2:2:end);
Sqx = zeros(4,size(elementsFEM,1)); Sqy = zeros(4,size(elementsFEM,1)); Sqc = zeros(4,size(elementsFEM,1));
for j = 1:size(elementsFEM,1)
    Sqx(1:4,j) = coordinatesFEM(elementsFEM(j,1:4),1);
    Sqy(1:4,j) = coordinatesFEM(elementsFEM(j,1:4),2);
    Sqc(1:4,j) = Uy(elementsFEM(j,1:4));
end
if size(elementsFEM,1)>2e4
    patch(Sqx,Sqy,Sqc,'facecolor','interp','edgecolor','none');
else
    patch(Sqx,Sqy,Sqc,'facecolor','interp');
end
view(2); axis tight; axis equal; colorbar; xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
set(gcf,'color','w'); colormap jet;
title('Click bad y-disp points and press -Enter- after done','fontweight','normal')
% figure; surf(v); colorbar; view(2)
% title('Displacement v','fontweight','normal')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Do you clear bad points by directly pointing y-disp bad points? (0-yes; 1-no)  \n')
prompt = 'Input here: ';
ClearBadInitialPointsOrNot = input(prompt);

while ClearBadInitialPointsOrNot == 0
    
    [row1, col1] = ginput; row = round(col1); col = round(row1);
    
    % Make sure it's within the image domain
    BadptCol=[row]; BadptRow=[col];
    [tempindex1] = find(BadptRow<1+max(coordinatesFEM(:,1)));  [tempindex2] = find(BadptCol<1+max(coordinatesFEM(:,2)));
    [tempindex3] = find(BadptRow>min(coordinatesFEM(:,1)));  [tempindex4] = find(BadptCol>min(coordinatesFEM(:,2)));
    tempindex1 = reshape(tempindex1,length(tempindex1),1);
    tempindex2 = reshape(tempindex2,length(tempindex2),1);
    tempindex3 = reshape(tempindex3,length(tempindex3),1);
    tempindex4 = reshape(tempindex4,length(tempindex4),1);
    tempindex = unique(intersect(tempindex4,intersect(tempindex3,intersect(tempindex1,tempindex2))));
    
    BadptRow = BadptRow(tempindex); BadptCol = BadptCol(tempindex);
    row=BadptRow; col=BadptCol;
    
    
    % --- Find elements include clicked bad points ---
    elementFEMCenterCoordx = 0.25* (coordinatesFEM( elementsFEM(:,1)  , 1) + coordinatesFEM( elementsFEM(:,2)  , 1) + ...
        coordinatesFEM( elementsFEM(:,3)  , 1) + coordinatesFEM( elementsFEM(:,4)  , 1));
    elementFEMCenterCoordy = 0.25* (coordinatesFEM( elementsFEM(:,1)  , 2) + coordinatesFEM( elementsFEM(:,2)  , 2) + ...
        coordinatesFEM( elementsFEM(:,3)  , 2) + coordinatesFEM( elementsFEM(:,4)  , 2));
    
    BadptCoord = [];
    for tempk = 1:length(row) % iterate for each clicking points
        
        DistElementAll = sqrt( ( elementFEMCenterCoordx-row(tempk) ).^2 + (elementFEMCenterCoordy-col(tempk)).^2 );
        DistCoordAll = sqrt( ( coordinatesFEM(:,1)-row(tempk) ).^2 + (coordinatesFEM(:,2)-col(tempk)).^2 );
        [row1,~] = find(DistElementAll < 1.42*mean(winstepsize));
        [row2,~] = find(DistCoordAll < 1.42*mean(winstepsize));
        temp1 = elementsFEM(row1,:); temp2 = row2(:);
        
        BadptCoord = [BadptCoord; unique([temp1(:);temp2(:)])]; BadptCoord=BadptCoord(:);
        
    end
    % Set unique BadptCoord
    BadptCoord = unique(BadptCoord); BadptCoord = setdiff(BadptCoord,[0]);
    
    U(2*BadptCoord-1) = NaN; U(2*BadptCoord) = NaN;
    F(4*BadptCoord-3) = NaN; F(4*BadptCoord-2) = NaN; F(4*BadptCoord-1) = NaN; F(4*BadptCoord) = NaN;
    
    %%%%%% Fill nans %%%%%%
    nanindex = find(isnan(U(1:2:end))==1); notnanindex = setdiff([1:1:size(coordinatesFEM,1)],nanindex);
    nanindexF = find(isnan(F(1:4:end))==1); notnanindexF = setdiff([1:1:size(coordinatesFEM,1)],nanindexF);
    
    % dilatedI = ( imgaussfilt(double(Df.ImgRefMask),0.5) );
    %  dilatedI = logical( dilatedI > 0.01); % figure, imshow(dilatedI)
    dilatedI = logical(DICpara.ImgRefMask);
    cc = bwconncomp(dilatedI,8);
    indPxAll = sub2ind( DICpara.ImgSize, round(coordinatesFEM(:,1)), round(coordinatesFEM(:,2)));
    indPxNotNanAll = sub2ind( DICpara.ImgSize, round(coordinatesFEM(notnanindex,1)), round(coordinatesFEM(notnanindex,2)) );
    stats = regionprops(cc,'Area','PixelList');
    
    for tempi = 1:length(stats)
        %%%%% Find those nodes %%%%%
        indPxtempi = sub2ind( DICpara.ImgSize, stats(tempi).PixelList(:,2), stats(tempi).PixelList(:,1) );
        Lia = ismember(indPxAll,indPxtempi); [LiaList,~] = find(Lia==1);
        Lib = ismember(indPxNotNanAll,indPxtempi); [LibList,~] = find(Lib==1);
        
        %%%%% RBF (Radial basis function) works better than "scatteredInterpolant" %%%%%
        % ------ Disp u ------
        op1 = rbfcreate( [coordinatesFEM(notnanindex(LibList),1:2)]',[U(2*notnanindex(LibList)-1)]','RBFFunction', 'thinplate'); rbfcheck(op1);
        fi1 = rbfinterp( [coordinatesFEM(LiaList,1:2)]', op1);
        U(2*LiaList-1) = fi1(:);
        % figure, plot3(coordinatesFEM(notnanindex(LibList),1),coordinatesFEM(notnanindex(LibList),2),U(2*notnanindex(LibList)-1),'.')
        % hold on; plot3(coordinatesFEM(LiaList,1),coordinatesFEM(LiaList,2),fi1,'.')
        
        % ------ Disp v ------
        op1 = rbfcreate( [coordinatesFEM(notnanindex(LibList),1:2)]',[U(2*notnanindex(LibList) )]','RBFFunction', 'thinplate'); rbfcheck(op1);
        fi1 = rbfinterp( [coordinatesFEM(LiaList,1:2)]', op1);
        U(2*LiaList ) = fi1(:);
        % ------ F11 ------
        op1 = rbfcreate( [coordinatesFEM(notnanindex(LibList),1:2)]',[F(4*notnanindex(LibList)-3)]','RBFFunction', 'thinplate'); rbfcheck(op1);
        fi1 = rbfinterp( [coordinatesFEM(LiaList,1:2)]', op1);
        F(4*LiaList-3) = fi1(:);
        % ------ F21 ------
        op1 = rbfcreate( [coordinatesFEM(notnanindex(LibList),1:2)]',[F(4*notnanindex(LibList)-2)]','RBFFunction', 'thinplate'); rbfcheck(op1);
        fi1 = rbfinterp( [coordinatesFEM(LiaList,1:2)]', op1);
        F(4*LiaList-2) = fi1(:);
        % ------ F12 ------
        op1 = rbfcreate( [coordinatesFEM(notnanindex(LibList),1:2)]',[F(4*notnanindex(LibList)-1)]','RBFFunction', 'thinplate'); rbfcheck(op1);
        fi1 = rbfinterp( [coordinatesFEM(LiaList,1:2)]', op1);
        F(4*LiaList-1) = fi1(:);
        % ------ F22 ------
        op1 = rbfcreate( [coordinatesFEM(notnanindex(LibList),1:2)]',[F(4*notnanindex(LibList)-0)]','RBFFunction', 'thinplate'); rbfcheck(op1);
        fi1 = rbfinterp( [coordinatesFEM(LiaList,1:2)]', op1);
        F(4*LiaList) = fi1(:);
    end
    
    %%%%% Replace these old codes by applying RBF interpolations %%%%%
    %     % ------ inpaint nans using gridfit ------
    %     Coordxnodes = unique(coordinatesFEM(:,1)); Coordynodes = unique(coordinatesFEM(:,2));
    %     nanindex = find(isnan(U(1:2:end))==1); notnanindex = setdiff([1:1:size(coordinatesFEM,1)],nanindex);
    %     
    %     
    %     %[CoordxnodesGrid,CoordynodesGrid] = ndgrid(Coordxnodes,Coordynodes);
    %     Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), U(2*notnanindex-1),'natural','nearest');
    %     U1 = Ftemp(coordinatesFEM(:,1),coordinatesFEM(:,2));
    %     %[u1temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %u1temp = u1temp';
    %     Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), U(2*notnanindex),'natural','nearest');
    %     V1 = Ftemp(coordinatesFEM(:,1),coordinatesFEM(:,2));
    %     %[v1temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %v1temp = v1temp';
    %     Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-3),'natural','nearest');
    %     F11 = Ftemp(coordinatesFEM(:,1),coordinatesFEM(:,2));
    %     %[F11temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %F11temp = F11temp';
    %     Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-2),'natural','nearest');
    %     F21 = Ftemp(coordinatesFEM(:,1),coordinatesFEM(:,2));
    %     %[F21temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %F21temp = F21temp';
    %     Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-1),'natural','nearest');
    %     F12 = Ftemp(coordinatesFEM(:,1),coordinatesFEM(:,2));
    %     %[F12temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %F12temp = F12temp';
    %     Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-0),'natural','nearest');
    %     F22 = Ftemp(coordinatesFEM(:,1),coordinatesFEM(:,2));
    %     %[F22temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %F22temp = F22temp';
    %     
    %     U = [U1(:),V1(:)]'; U = U(:);
    %     F = [F11(:),F21(:),F12(:),F22(:)]'; F = F(:);
    
    
    close all;  Plotdisp_show(U,coordinatesFEM,elementsFEM);
    Uy = U(2:2:end);
    Sqx = zeros(4,size(elementsFEM,1)); Sqy = zeros(4,size(elementsFEM,1)); Sqc = zeros(4,size(elementsFEM,1));
    for j = 1:size(elementsFEM,1)
        Sqx(1:4,j) = coordinatesFEM(elementsFEM(j,1:4),1);
        Sqy(1:4,j) = coordinatesFEM(elementsFEM(j,1:4),2);
        Sqc(1:4,j) = Uy(elementsFEM(j,1:4));
    end
    if size(elementsFEM,1)>2e4
        patch(Sqx,Sqy,Sqc,'facecolor','interp','edgecolor','none');
    else
        patch(Sqx,Sqy,Sqc,'facecolor','interp');
    end
    view(2); axis tight; axis equal; colorbar; xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
    set(gcf,'color','w'); colormap jet;
    title('Click bad dispy-v points and press -Enter- after done','fontweight','normal')
    
    prompt = 'Do you point out more y-disp bad points? (0-yes; 1-no) Input: ';
    ClearBadInitialPointsOrNot = input(prompt);
    
end




%% Remove bad F values
%%%%%% Fill nans %%%%%%
nanindex = find(isnan(U(1:2:end))==1); notnanindex = setdiff([1:1:size(coordinatesFEM,1)],nanindex);
nanindexF = find(isnan(F(1:4:end))==1); notnanindexF = setdiff([1:1:size(coordinatesFEM,1)],nanindexF);

if ~isempty(nanindex) || ~isempty(nanindexF)
    
Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1),coordinatesFEM(notnanindex,2),U(2*notnanindex-1),'nearest','nearest');
U1 = Ftemp(coordinatesFEM(:,1),coordinatesFEM(:,2));
Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1),coordinatesFEM(notnanindex,2),U(2*notnanindex),'nearest','nearest');
U2 = Ftemp(coordinatesFEM(:,1),coordinatesFEM(:,2));
Ftemp = scatteredInterpolant(coordinatesFEM(notnanindexF,1),coordinatesFEM(notnanindexF,2),F(4*notnanindexF-3),'nearest','nearest');
F11 = Ftemp(coordinatesFEM(:,1),coordinatesFEM(:,2));
Ftemp = scatteredInterpolant(coordinatesFEM(notnanindexF,1),coordinatesFEM(notnanindexF,2),F(4*notnanindexF-2),'nearest','nearest');
F21 = Ftemp(coordinatesFEM(:,1),coordinatesFEM(:,2));
Ftemp = scatteredInterpolant(coordinatesFEM(notnanindexF,1),coordinatesFEM(notnanindexF,2),F(4*notnanindexF-1),'nearest','nearest');
F12 = Ftemp(coordinatesFEM(:,1),coordinatesFEM(:,2));
Ftemp = scatteredInterpolant(coordinatesFEM(notnanindexF,1),coordinatesFEM(notnanindexF,2),F(4*notnanindexF-0),'nearest','nearest');
F22 = Ftemp(coordinatesFEM(:,1),coordinatesFEM(:,2));

U = [U1(:),U2(:)]'; U = U(:);
F = [F11(:),F21(:),F12(:),F22(:)]'; F = F(:);

end

% ------ inpaint nans using gridfit ------
% minCoordStep = min( [DICmesh.elementMinSize] );
% xList = [ceil(min(coordinatesFEM(:,1))) : minCoordStep : floor(max(coordinatesFEM(:,1)))]';
% yList = [ceil(min(coordinatesFEM(:,2))) : minCoordStep : floor(max(coordinatesFEM(:,2)))]';
% 
% [xGrid,yGrid] = ndgrid(xList,yList);
% smoothness = 1e-4;
% tempF11 = regularizeNd([coordinatesFEM(:,1), coordinatesFEM(:,2)],F(1:4:end),{xList,yList},smoothness);
% tempF21 = regularizeNd([coordinatesFEM(:,1), coordinatesFEM(:,2)],F(2:4:end),{xList,yList},smoothness);
% tempF12 = regularizeNd([coordinatesFEM(:,1), coordinatesFEM(:,2)],F(3:4:end),{xList,yList},smoothness);
% tempF22 = regularizeNd([coordinatesFEM(:,1), coordinatesFEM(:,2)],F(4:4:end),{xList,yList},smoothness);
% % tempU1 = regularizeNd([coordinatesFEM(:,1), coordinatesFEM(:,2)],U(1:2:end),{xList,yList},smoothness/10);
% % tempU2 = regularizeNd([coordinatesFEM(:,1), coordinatesFEM(:,2)],U(1:2:end),{xList,yList},smoothness/10);
% 
% F_F11 = scatteredInterpolant(xGrid(:),yGrid(:),tempF11(:),'natural','nearest');
% F11 = F_F11(coordinatesFEM(:,1),coordinatesFEM(:,2));
% F_F21 = scatteredInterpolant(xGrid(:),yGrid(:),tempF21(:),'natural','nearest');
% F21 = F_F21(coordinatesFEM(:,1),coordinatesFEM(:,2));
% F_F12 = scatteredInterpolant(xGrid(:),yGrid(:),tempF12(:),'natural','nearest');
% F12 = F_F12(coordinatesFEM(:,1),coordinatesFEM(:,2));
% F_F22 = scatteredInterpolant(xGrid(:),yGrid(:),tempF22(:),'natural','nearest');
% F22 = F_F22(coordinatesFEM(:,1),coordinatesFEM(:,2));
% 
% F = [F11(:),F21(:),F12(:),F22(:)]'; F = F(:);

 

end




%% ========================================================================
function [medianU, normFluct] = funMedRemoveOutliers(u,epsilon)

inpaint_opt = 3;

nSize = 2*[1 1];
skipIdx = ceil(prod(nSize)/2);
padOption = 'symmetric';

u =  inpaint_nans(double(u),inpaint_opt);

medianU = medFilt2(u,nSize,padOption,skipIdx);
fluct = u - medianU;
medianRes = medFilt2(abs(fluct),nSize,padOption,skipIdx);
normFluct = abs(fluct./(medianRes + epsilon));


end

%% ========================================================================
function Vr = medFilt2(V0,nSize, padoption, skipIdx)
% fast median filter for 2D data with extra options.

if nargin < 4, skipIdx = 0; end
if nargin < 3, padoption = 'symmetric'; end
if nargin < 2, nSize = [2 2]; end

nLength = prod(nSize);
if mod(nLength,2) == 1, padSize = floor(nSize/2);
elseif mod(nLength,2) == 0, padSize = [nSize(1)/2-1,nSize(2)/2];
end

if strcmpi(padoption,'none')
    V = V0;
else
    V = (padarray(V0,padSize(1)*[1,1],padoption,'pre'));
    V = (padarray(V,padSize(2)*[1,1],padoption,'post'));
end

S = size(V);
nLength = prod(nSize)-sum(skipIdx>1);
Vn = single(zeros(S(1)-(nSize(1)-1),S(2)-(nSize(2)-1),nLength));  % all the neighbor

%%
% build the neighborhood

i = cell(1,nSize(1)); j = cell(1,nSize(2));
for m = 1:nSize(1), i{m} = m:(S(1)-(nSize(1)-m)); end
for m = 1:nSize(2), j{m} = m:(S(2)-(nSize(2)-m)); end

p = 1;
for m = 1:nSize(1)
    for n = 1:nSize(2)
        if p ~= skipIdx || skipIdx == 0
            Vn(:,:,p) = V(i{m},j{n});
        end
        p = p + 1;
    end
end

if skipIdx ~= 0, Vn(:,:,skipIdx) = []; end
% perform the processing
Vn = sort(Vn,3);

if mod(nLength,2) == 1 % if odd get the middle element
    Vr = Vn(:,:,ceil(nLength/2));
else % if even get the mean of the two middle elements
    Vr = mean(cat(3,Vn(:,:,nLength/2),Vn(:,:,nLength/2+1)),4);
end

end

%% ========================================================================
function [cc, ccMask] = ...
    removeBadCorrelations(cc,ccThreshold,sizeChange,mSize)

if sizeChange == 1
    %recompute threshold, only use pce & ppe since these give the best
    %results emprically.
    for ii = 1:2
        
        [qf_para{ii},single_distro] = bimodal_gauss_fit(cc.qfactors(:,ii));
        
        if single_distro == 0%(qf_para{ii}(2) + 2*qf_para{ii}(4)) < (qf_para{ii}(3) - 2*qf_para{ii}(5))
            cc.q_thresh{ii} = qf_para{ii}(3) - ccThreshold*qf_para{ii}(5);
        elseif single_distro == 1
            cc.q_thresh{ii} = qf_para{ii}(3) - ccThreshold*qf_para{ii}(5);
        else
            cc.q_thresh{ii} = qf_para{ii}(3) - ccThreshold*qf_para{ii}(5);
        end
    end
    q_trim = [cc.q_thresh{1};cc.q_thresh{2}];
else
    q_trim = [cc.q_thresh{1};cc.q_thresh{2}];
end

%NaN the qfactor values that are below the threshold
temp = bsxfun(@le,cc.qfactors(:,1:2),q_trim');
qfactors_accept = cc.qfactors(:,1:2);
qfactors_accept(temp) = NaN;

for ii = 1:2
    cc.qfactors_accept{ii} = reshape(double(qfactors_accept(:,ii)),mSize);
end

ccMask = ones(size(qfactors_accept)) + ...
    zeros(size(qfactors_accept)).*qfactors_accept;

end



%% ========================================================================
function [paramEsts,single_distro] = bimodal_gauss_fit(x)
%This function takes a dataset and fits a bimodal Gaussian distro to it.

x = sort(x);

%set function for bimodal Gaussian
pdf_normmixture = @(x,p,mu1,mu2,sigma1,sigma2) ...
    p*normpdf(x,mu1,sigma1) + (1-p)*normpdf(x,mu2,sigma2);
pdf_single = @(x,mu1,sigma1) ...
    normpdf(x,mu1,sigma1);

%starting params, biased mixture toward "good" values,
%centered at quartiles, equal std dev.
pStart = 0.25;
muStart = quantile(x,[.10 .75]);
sigmaStart(1) = sqrt(var(x(1:round(length(x)/5))));
%- 0.25*diff(quantile(x,[0.01 0.25])).^2);
sigmaStart(2) = sqrt(var(x(ceil(length(x)/10):ceil(3*length(x)/4))));
%... - 0.25*diff(quantile(x,[0.25 0.75])).^2);%1:round(length(x)/2)
start = [pStart muStart sigmaStart];

%set lower and upper bounds
lb = [0 -inf -inf 0.00001 0.00001];
ub = [1 inf inf inf inf];

%do the parameter estimation
options = statset('MaxIter',1800, 'MaxFunEvals',3600);
% options.FunValCheck = 'off';
try
    single_distro = 0;
    paramEsts = mle(x, 'pdf',pdf_normmixture, 'start',start, ...
        'lower',lb, 'upper',ub, 'options',options);%,'optimfun','fmincon'
    
    if paramEsts(2)-paramEsts(4) >= paramEsts(3)+paramEsts(5) || ...
            paramEsts(2)+paramEsts(4) <= paramEsts(3)-paramEsts(5)
        
        single_distro = 1;
        %     disp('Parameters estimated for single peak Gaussian')
        paramEsts = mle(x,'options',options);%,'optimfun','fmincon'
        paramEsts = [0.5,paramEsts(1),paramEsts(1),paramEsts(2),...
            paramEsts(2)];
        
    end
    
catch
    single_distro = 1;
    %     disp('Parameters estimated for single peak Gaussian')
    paramEsts = mle(x,'options',options);%,'optimfun','fmincon'
    paramEsts = [0.5,paramEsts(1),paramEsts(1),paramEsts(2),...
        paramEsts(2)];
end

% % %show the result
% % figure
% % % [~, bins] =
% % histogram(x,100);
% % % bins = -2.5:.5:7.5;
% % % h = bar(bins,histc(x,bins)/(length(x)*0.5),'histc');
% % % histogram(x,100)
% % % h.FaceColor = [0.9 0.9 0.9];
% % xgrid = linspace(1.1*min(x),1.1*max(x),200);
% % pdfgrid = pdf_normmixture(xgrid,paramEsts(1),paramEsts(2),paramEsts(3),...
% %     paramEsts(4),paramEsts(5));
% % hold on
% % plot((paramEsts(3) - 2*paramEsts(5)),pdfgrid,'or')
% % plot((paramEsts(2) + 2*paramEsts(4)),pdfgrid,'*r')
% % plot(xgrid,pdfgrid,'-b')
% % hold off
% % xlabel('x')
% % ylabel('Probability Density')

end




