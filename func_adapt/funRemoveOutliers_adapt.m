function [U,F] = funRemoveOutliers_adapt(coordinatesFEM,elementsFEM,U,F,Thr0,varargin)
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


switch nargin
    case 7
        BadptRow = varargin{1}; BadptCol = varargin{2};
    otherwise
        BadptRow = []; BadptCol = [];
end

%% Setting upper and lower bounds
close all; %Plotdisp_show(U,coordinatesFEM,elementsFEM);
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
% figure; surf(v); colorbar; view(2)
% title('Displacement v','fontweight','normal')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ==============================================
% Manually remove bad points.
% if Thr0 ~= 100
    
fprintf('Do you clear bad points by setting upper/lower bounds once more? (0-yes; 1-no)  \n')
prompt = 'Input here: ';
ClearBadInitialPointsOrNot = input(prompt);

while ClearBadInitialPointsOrNot == 0
    
    prompt = 'What is your upper bound for x-displacement? Input: ';
    upperbound = input(prompt);
    [row1,~] = find(U(1:2:end)>upperbound);
    prompt = 'What is your lower bound for x-displacement? Input: ';
    lowerbound = input(prompt);
    [row2,~] = find(U(1:2:end)<lowerbound);
    prompt = 'What is your upper bound for y-displacement? Input: ';
    upperbound = input(prompt);
    [row3,~] = find(U(2:2:end)>upperbound);
    prompt = 'What is your lower bound for y-displacement? Input: ';
    lowerbound = input(prompt);
    [row4,~] = find(U(2:2:end)<lowerbound);
     
    row = unique([row1; row2; row3; row4]);  
    
    U(2*row-1) = NaN; U(2*row) = NaN; 
    F(4*row-3) = NaN; F(4*row-2) = NaN; F(4*row-1) = NaN; F(4*row) = NaN;

    % ------ inpaint nans using gridfit ------
    Coordxnodes = unique(coordinatesFEM(:,1)); Coordynodes = unique(coordinatesFEM(:,2));
    nanindex = find(isnan(U(1:2:end))==1); notnanindex = setdiff([1:1:size(coordinatesFEM,1)],nanindex);
    
    %[u1temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), U(2*notnanindex-1),Coordxnodes,Coordynodes,'regularizer','springs'); u1temp = u1temp';
    %[v1temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), U(2*notnanindex),Coordxnodes,Coordynodes,'regularizer','springs'); v1temp = v1temp';
    %[F11temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-3),Coordxnodes,Coordynodes,'regularizer','springs'); F11temp = F11temp';
    %[F21temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-2),Coordxnodes,Coordynodes,'regularizer','springs'); F21temp = F21temp';
    %[F12temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-1),Coordxnodes,Coordynodes,'regularizer','springs'); F12temp = F12temp';
    %[F22temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-0),Coordxnodes,Coordynodes,'regularizer','springs'); F22temp = F22temp';

    [CoordxnodesGrid,CoordynodesGrid] = ndgrid(Coordxnodes,Coordynodes);
    Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), U(2*notnanindex-1));
    [u1temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %u1temp = u1temp';
    Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), U(2*notnanindex));
    [v1temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %v1temp = v1temp';
    Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-3));
    [F11temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %F11temp = F11temp';
    Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-2));
    [F21temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %F21temp = F21temp';
    Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-1));
    [F12temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %F12temp = F12temp';
    Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-0));
    [F22temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %F22temp = F22temp';
     
    for tempi = 1:size(coordinatesFEM,1)
        [row1,col1] = find(Coordxnodes==coordinatesFEM(tempi,1));
        [row2,col2] = find(Coordynodes==coordinatesFEM(tempi,2));
        U(2*tempi-1) = u1temp(row1,row2);
        U(2*tempi)   = v1temp(row1,row2);
        F(4*tempi-3) = F11temp(row1,row2); F(4*tempi-2) = F21temp(row1,row2);
        F(4*tempi-1) = F12temp(row1,row2); F(4*tempi) = F22temp(row1,row2);
    end

    close all;  
    Ux = U(1:2:end); figure,
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
    
    
    Uy = U(2:2:end); figure,
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
    
     prompt = 'Do you want to reset upper/lower bounds? (0-yes; 1-no) Input: ';
    ClearBadInitialPointsOrNot = input(prompt);
    
end

    







%%
%coordinatesFEM = coordinatesFEMIter; elementsFEM = elementsFEMIter;
%U = UIter; F = FIter;
close all; %Plotdisp_show(U,coordinatesFEM,elementsFEM);
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
 

%% =========
fprintf('Do you clear bad points by directly pointing x-disp bad points? (0-yes; 1-no)  \n')
prompt = 'Input here: ';
ClearBadInitialPointsOrNot = input(prompt);

while ClearBadInitialPointsOrNot == 0
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Have a look at integer search
    % --------------------------------------
%     close all; %Plotdisp_show(U,coordinatesFEM,elementsFEM);
%     Ux = U(1:2:end);
%     Sqx = zeros(4,size(elementsFEM,1)); Sqy = zeros(4,size(elementsFEM,1)); Sqc = zeros(4,size(elementsFEM,1));
%     for j = 1:size(elementsFEM,1)
%         Sqx(1:4,j) = coordinatesFEM(elementsFEM(j,1:4),1);
%         Sqy(1:4,j) = coordinatesFEM(elementsFEM(j,1:4),2);
%         Sqc(1:4,j) = Ux(elementsFEM(j,1:4));
%     end
%     patch(Sqx,Sqy,Sqc,'facecolor','interp','edgecolor','none');
%     view(2); axis tight; axis equal; colorbar; xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
%     set(gcf,'color','w'); colormap jet;
%     title('Click bad dispx-u points and press -Enter- after done','fontweight','normal')
    % figure; surf(v); colorbar; view(2)
    % title('Displacement v','fontweight','normal')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [row1, col1] = ginput; row = floor(col1); col = floor(row1); 
    %row13 = row; col13 = col+2;     row22 = row-1; col22 = col+1;
    %row23 = row; col23 = col+1;     row24 = row+1; col24 = col+1;
    %row31 = row-2; col31 = col;     row32 = row-1; col32 = col; 
    %row34 = row+1; col34 = col;     row35 = row+2; col35 = col;
    %row42 = row-1; col42 = col-1;   row43 = row; col43 = col-1;
    %row44 = row+1; col44 = col-1;   row53 = row; col53 = col-2;
    %row = [row;row13;row22;row23;row24;row31;row32;row34;row35;row42;row43;row44;row53];
    %col = [col;col13;col22;col23;col24;col31;col32;col34;col35;col42;col43;col44;col53];
     
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
     
     % --- Find which elements ---
     BadptCoord = [];
    for tempk = 1:length(row) % iterate for each clicking points
        % Compute distance
        DistCoordAll = (coordinatesFEM(:,2)-col(tempk)).^2 + (coordinatesFEM(:,1)-row(tempk)).^2; 
        [DistCoordAllSorted, DistCoordAllSortedInd] = sort(sqrt(DistCoordAll));
        ismember1 = []; ismember2 = []; ismember3 = []; 
        for templ = 1:size(elementsFEM,1)
            ismember1(templ) = ismember(DistCoordAllSortedInd(1),elementsFEM(templ,1:8));
            ismember2(templ) = ismember(DistCoordAllSortedInd(2),elementsFEM(templ,1:8));
            ismember3(templ) = ismember(DistCoordAllSortedInd(3),elementsFEM(templ,1:8));
        end
        [~,coll] = find(ismember1+ismember2+ismember3 == 3);
        BadptCoord = [BadptCoord;elementsFEM(coll,:)'];
    end
    % Set unique BadptCoord
    BadptCoord = unique(BadptCoord); BadptCoord = setdiff(BadptCoord,[0]);
    
    U(2*BadptCoord-1) = NaN; U(2*BadptCoord) = NaN; 
    F(4*BadptCoord-3) = NaN; F(4*BadptCoord-2) = NaN; F(4*BadptCoord-1) = NaN; F(4*BadptCoord) = NaN;

    % ------ inpaint nans using gridfit ------
    Coordxnodes = unique(coordinatesFEM(:,1)); Coordynodes = unique(coordinatesFEM(:,2));
    nanindex = find(isnan(U(1:2:end))==1); notnanindex = setdiff([1:1:size(coordinatesFEM,1)],nanindex);
    
    %[u1temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), U(2*notnanindex-1),Coordxnodes,Coordynodes,'regularizer','springs'); u1temp = u1temp';
    %[v1temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), U(2*notnanindex),Coordxnodes,Coordynodes,'regularizer','springs'); v1temp = v1temp';
    %[F11temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-3),Coordxnodes,Coordynodes,'regularizer','springs'); F11temp = F11temp';
    %[F21temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-2),Coordxnodes,Coordynodes,'regularizer','springs'); F21temp = F21temp';
    %[F12temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-1),Coordxnodes,Coordynodes,'regularizer','springs'); F12temp = F12temp';
    %[F22temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-0),Coordxnodes,Coordynodes,'regularizer','springs'); F22temp = F22temp';

    [CoordxnodesGrid,CoordynodesGrid] = ndgrid(Coordxnodes,Coordynodes);
    Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), U(2*notnanindex-1));
    [u1temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %u1temp = u1temp';
    Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), U(2*notnanindex));
    [v1temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %v1temp = v1temp';
    Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-3));
    [F11temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %F11temp = F11temp';
    Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-2));
    [F21temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %F21temp = F21temp';
    Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-1));
    [F12temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %F12temp = F12temp';
    Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-0));
    [F22temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %F22temp = F22temp';
     
    for tempi = 1:size(coordinatesFEM,1)
        [row1,col1] = find(Coordxnodes==coordinatesFEM(tempi,1));
        [row2,col2] = find(Coordynodes==coordinatesFEM(tempi,2));
        U(2*tempi-1) = u1temp(row1,row2);
        U(2*tempi)   = v1temp(row1,row2);
        F(4*tempi-3) = F11temp(row1,row2); F(4*tempi-2) = F21temp(row1,row2);
        F(4*tempi-1) = F12temp(row1,row2); F(4*tempi) = F22temp(row1,row2);
    end

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




%% =========
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Have a look at integer search
% --------------------------------------
close all; %Plotdisp_show(U,coordinatesFEM,elementsFEM);
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
% figure; surf(v); colorbar; view(2)
% title('Displacement v','fontweight','normal')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Do you clear bad points by directly pointing y-disp bad points? (0-yes; 1-no)  \n')
prompt = 'Input here: ';
ClearBadInitialPointsOrNot = input(prompt);

while ClearBadInitialPointsOrNot == 0
     
    [row1, col1] = ginput; row = floor(col1); col = floor(row1); 
    %row13 = row; col13 = col+2;     row22 = row-1; col22 = col+1;
    %row23 = row; col23 = col+1;     row24 = row+1; col24 = col+1;
    %row31 = row-2; col31 = col;     row32 = row-1; col32 = col; 
    %row34 = row+1; col34 = col;     row35 = row+2; col35 = col;
    %row42 = row-1; col42 = col-1;   row43 = row; col43 = col-1;
    %row44 = row+1; col44 = col-1;   row53 = row; col53 = col-2;
    %row = [row;row13;row22;row23;row24;row31;row32;row34;row35;row42;row43;row44;row53];
    %col = [col;col13;col22;col23;col24;col31;col32;col34;col35;col42;col43;col44;col53];
     
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
     
     % --- Find which elements ---
     BadptCoord = [];
    for tempk = 1:length(row) % iterate for each clicking points
        % Compute distance
        DistCoordAll = (coordinatesFEM(:,2)-col(tempk)).^2 + (coordinatesFEM(:,1)-row(tempk)).^2; 
        [DistCoordAllSorted, DistCoordAllSortedInd] = sort(sqrt(DistCoordAll));
        ismember1 = []; ismember2 = []; ismember3 = []; 
        for templ = 1:size(elementsFEM,1)
            ismember1(templ) = ismember(DistCoordAllSortedInd(1),elementsFEM(templ,1:8));
            ismember2(templ) = ismember(DistCoordAllSortedInd(2),elementsFEM(templ,1:8));
            ismember3(templ) = ismember(DistCoordAllSortedInd(3),elementsFEM(templ,1:8));
        end
        [~,coll] = find(ismember1+ismember2+ismember3 == 3);
        BadptCoord = [BadptCoord;elementsFEM(coll,:)'];
    end
    % Set unique BadptCoord
    BadptCoord = unique(BadptCoord); BadptCoord = setdiff(BadptCoord,[0]);
    
    U(2*BadptCoord-1) = NaN; U(2*BadptCoord) = NaN; 
    F(4*BadptCoord-3) = NaN; F(4*BadptCoord-2) = NaN; F(4*BadptCoord-1) = NaN; F(4*BadptCoord) = NaN;

    % ------ inpaint nans using gridfit ------
    Coordxnodes = unique(coordinatesFEM(:,1)); Coordynodes = unique(coordinatesFEM(:,2));
    nanindex = find(isnan(U(1:2:end))==1); notnanindex = setdiff([1:1:size(coordinatesFEM,1)],nanindex);
%     [u1temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), U(2*notnanindex-1),Coordxnodes,Coordynodes,'regularizer','springs'); u1temp = u1temp';
%     [v1temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), U(2*notnanindex),Coordxnodes,Coordynodes,'regularizer','springs'); v1temp = v1temp';
%     [F11temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-3),Coordxnodes,Coordynodes,'regularizer','springs'); F11temp = F11temp';
%     [F21temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-2),Coordxnodes,Coordynodes,'regularizer','springs'); F21temp = F21temp';
%     [F12temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-1),Coordxnodes,Coordynodes,'regularizer','springs'); F12temp = F12temp';
%     [F22temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-0),Coordxnodes,Coordynodes,'regularizer','springs'); F22temp = F22temp';

    [CoordxnodesGrid,CoordynodesGrid] = ndgrid(Coordxnodes,Coordynodes);
    Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), U(2*notnanindex-1));
    [u1temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %u1temp = u1temp';
    Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), U(2*notnanindex));
    [v1temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %v1temp = v1temp';
    Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-3));
    [F11temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %F11temp = F11temp';
    Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-2));
    [F21temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %F21temp = F21temp';
    Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-1));
    [F12temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %F12temp = F12temp';
    Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-0));
    [F22temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %F22temp = F22temp';
     
    
    for tempi = 1:size(coordinatesFEM,1)
        [row1,col1] = find(Coordxnodes==coordinatesFEM(tempi,1));
        [row2,col2] = find(Coordynodes==coordinatesFEM(tempi,2));
        U(2*tempi-1) = u1temp(row1,row2);
        U(2*tempi)   = v1temp(row1,row2);
        F(4*tempi-3) = F11temp(row1,row2); F(4*tempi-2) = F21temp(row1,row2);
        F(4*tempi-1) = F12temp(row1,row2); F(4*tempi) = F22temp(row1,row2);
    end

    close all; %Plotdisp_show(U,coordinatesFEM,elementsFEM);
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




%%


% %% ========= F11 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Have a look at integer search
% % --------------------------------------
% close all; %Plotdisp_show(U,coordinatesFEM,elementsFEM);
% Ftemp = F(1:4:end);
% Sqx = zeros(4,size(elementsFEM,1)); Sqy = zeros(4,size(elementsFEM,1)); Sqc = zeros(4,size(elementsFEM,1));
% for j = 1:size(elementsFEM,1)
%     Sqx(1:4,j) = coordinatesFEM(elementsFEM(j,1:4),1);
%     Sqy(1:4,j) = coordinatesFEM(elementsFEM(j,1:4),2);
%     Sqc(1:4,j) = Ftemp(elementsFEM(j,1:4));
% end
% if size(elementsFEM,1)>2e4
%     patch(Sqx,Sqy,Sqc,'facecolor','interp','edgecolor','none'); 
% else
%     patch(Sqx,Sqy,Sqc,'facecolor','interp'); 
% end
% view(2); axis tight; axis equal; colorbar; xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
% set(gcf,'color','w'); colormap jet;
% title('Click bad F11 points and press -Enter- after done','fontweight','normal')
% % figure; surf(v); colorbar; view(2)
% % title('Displacement v','fontweight','normal')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fprintf('Do you clear bad points by directly pointing F11-bad points? (0-yes; 1-no)  \n')
% prompt = 'Input here: ';
% ClearBadInitialPointsOrNot = input(prompt);
% 
% while ClearBadInitialPointsOrNot == 0
%      
%     [row1, col1] = ginput; row = floor(col1); col = floor(row1); 
%     %row13 = row; col13 = col+2;     row22 = row-1; col22 = col+1;
%     %row23 = row; col23 = col+1;     row24 = row+1; col24 = col+1;
%     %row31 = row-2; col31 = col;     row32 = row-1; col32 = col; 
%     %row34 = row+1; col34 = col;     row35 = row+2; col35 = col;
%     %row42 = row-1; col42 = col-1;   row43 = row; col43 = col-1;
%     %row44 = row+1; col44 = col-1;   row53 = row; col53 = col-2;
%     %row = [row;row13;row22;row23;row24;row31;row32;row34;row35;row42;row43;row44;row53];
%     %col = [col;col13;col22;col23;col24;col31;col32;col34;col35;col42;col43;col44;col53];
%      
%     BadptCol=[row]; BadptRow=[col]; 
%     [tempindex1] = find(BadptRow<1+max(coordinatesFEM(:,1)));  [tempindex2] = find(BadptCol<1+max(coordinatesFEM(:,2)));
%     [tempindex3] = find(BadptRow>min(coordinatesFEM(:,1)));  [tempindex4] = find(BadptCol>min(coordinatesFEM(:,2)));
%     tempindex1 = reshape(tempindex1,length(tempindex1),1);
%     tempindex2 = reshape(tempindex2,length(tempindex2),1);
%     tempindex3 = reshape(tempindex3,length(tempindex3),1);
%     tempindex4 = reshape(tempindex4,length(tempindex4),1);
%     tempindex = unique(intersect(tempindex4,intersect(tempindex3,intersect(tempindex1,tempindex2))));
%  
%     BadptRow = BadptRow(tempindex); BadptCol = BadptCol(tempindex); 
%     row=BadptRow; col=BadptCol;
%      
%      % --- Find which elements ---
%      BadptCoord = [];
%     for tempk = 1:length(row) % iterate for each clicking points
%         % Compute distance
%         DistCoordAll = (coordinatesFEM(:,2)-col(tempk)).^2 + (coordinatesFEM(:,1)-row(tempk)).^2; 
%         [DistCoordAllSorted, DistCoordAllSortedInd] = sort(sqrt(DistCoordAll));
%         ismember1 = []; ismember2 = []; ismember3 = []; 
%         for templ = 1:size(elementsFEM,1)
%             ismember1(templ) = ismember(DistCoordAllSortedInd(1),elementsFEM(templ,1:8));
%             ismember2(templ) = ismember(DistCoordAllSortedInd(2),elementsFEM(templ,1:8));
%             ismember3(templ) = ismember(DistCoordAllSortedInd(3),elementsFEM(templ,1:8));
%         end
%         [~,coll] = find(ismember1+ismember2+ismember3 == 3);
%         BadptCoord = [BadptCoord;elementsFEM(coll,:)'];
%     end
%     % Set unique BadptCoord
%     BadptCoord = unique(BadptCoord); BadptCoord = setdiff(BadptCoord,[0]);
%     
%     U(2*BadptCoord-1) = NaN; U(2*BadptCoord) = NaN; 
%     F(4*BadptCoord-3) = NaN; F(4*BadptCoord-2) = NaN; F(4*BadptCoord-1) = NaN; F(4*BadptCoord) = NaN;
% 
%     % ------ inpaint nans using gridfit ------
%     Coordxnodes = unique(coordinatesFEM(:,1)); Coordynodes = unique(coordinatesFEM(:,2));
%     nanindex = find(isnan(U(1:2:end))==1); notnanindex = setdiff([1:1:size(coordinatesFEM,1)],nanindex);
% %     [u1temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), U(2*notnanindex-1),Coordxnodes,Coordynodes,'regularizer','springs'); u1temp = u1temp';
% %     [v1temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), U(2*notnanindex),Coordxnodes,Coordynodes,'regularizer','springs'); v1temp = v1temp';
% %     [F11temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-3),Coordxnodes,Coordynodes,'regularizer','springs'); F11temp = F11temp';
% %     [F21temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-2),Coordxnodes,Coordynodes,'regularizer','springs'); F21temp = F21temp';
% %     [F12temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-1),Coordxnodes,Coordynodes,'regularizer','springs'); F12temp = F12temp';
% %     [F22temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-0),Coordxnodes,Coordynodes,'regularizer','springs'); F22temp = F22temp';
% 
%     [CoordxnodesGrid,CoordynodesGrid] = ndgrid(Coordxnodes,Coordynodes);
%     Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), U(2*notnanindex-1));
%     [u1temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %u1temp = u1temp';
%     Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), U(2*notnanindex));
%     [v1temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %v1temp = v1temp';
%     Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-3));
%     [F11temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %F11temp = F11temp';
%     Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-2));
%     [F21temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %F21temp = F21temp';
%     Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-1));
%     [F12temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %F12temp = F12temp';
%     Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-0));
%     [F22temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %F22temp = F22temp';
%      
%     
%     for tempi = 1:size(coordinatesFEM,1)
%         [row1,col1] = find(Coordxnodes==coordinatesFEM(tempi,1));
%         [row2,col2] = find(Coordynodes==coordinatesFEM(tempi,2));
%         U(2*tempi-1) = u1temp(row1,row2);
%         U(2*tempi)   = v1temp(row1,row2);
%         F(4*tempi-3) = F11temp(row1,row2); F(4*tempi-2) = F21temp(row1,row2);
%         F(4*tempi-1) = F12temp(row1,row2); F(4*tempi) = F22temp(row1,row2);
%     end
% 
%     close all; %Plotdisp_show(U,coordinatesFEM,elementsFEM);
%     Ftemp = F(1:4:end);
%     Sqx = zeros(4,size(elementsFEM,1)); Sqy = zeros(4,size(elementsFEM,1)); Sqc = zeros(4,size(elementsFEM,1));
%     for j = 1:size(elementsFEM,1)
%         Sqx(1:4,j) = coordinatesFEM(elementsFEM(j,1:4),1);
%         Sqy(1:4,j) = coordinatesFEM(elementsFEM(j,1:4),2);
%         Sqc(1:4,j) = Ftemp(elementsFEM(j,1:4));
%     end
%     if size(elementsFEM,1)>2e4
%         patch(Sqx,Sqy,Sqc,'facecolor','interp','edgecolor','none');
%     else
%         patch(Sqx,Sqy,Sqc,'facecolor','interp');
%     end
%     view(2); axis tight; axis equal; colorbar; xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
%     set(gcf,'color','w'); colormap jet;
%     title('Click bad F11 points and press -Enter- after done','fontweight','normal')
%  
%     prompt = 'Do you point out more F11-bad points? (0-yes; 1-no) Input: ';
%     ClearBadInitialPointsOrNot = input(prompt);
%     
% end
% 
% 
% 
% 
% %% ========= F21
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Have a look at integer search
% % --------------------------------------
% close all; %Plotdisp_show(U,coordinatesFEM,elementsFEM);
% Ftemp = F(2:4:end);
% Sqx = zeros(4,size(elementsFEM,1)); Sqy = zeros(4,size(elementsFEM,1)); Sqc = zeros(4,size(elementsFEM,1));
% for j = 1:size(elementsFEM,1)
%     Sqx(1:4,j) = coordinatesFEM(elementsFEM(j,1:4),1);
%     Sqy(1:4,j) = coordinatesFEM(elementsFEM(j,1:4),2);
%     Sqc(1:4,j) = Ftemp(elementsFEM(j,1:4));
% end
% if size(elementsFEM,1)>2e4
%     patch(Sqx,Sqy,Sqc,'facecolor','interp','edgecolor','none'); 
% else
%     patch(Sqx,Sqy,Sqc,'facecolor','interp'); 
% end
% view(2); axis tight; axis equal; colorbar; xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
% set(gcf,'color','w'); colormap jet;
% title('Click bad F21 points and press -Enter- after done','fontweight','normal')
% % figure; surf(v); colorbar; view(2)
% % title('Displacement v','fontweight','normal')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fprintf('Do you clear bad points by directly pointing F21-bad points? (0-yes; 1-no)  \n')
% prompt = 'Input here: ';
% ClearBadInitialPointsOrNot = input(prompt);
% 
% while ClearBadInitialPointsOrNot == 0
%      
%     [row1, col1] = ginput; row = floor(col1); col = floor(row1); 
%     %row13 = row; col13 = col+2;     row22 = row-1; col22 = col+1;
%     %row23 = row; col23 = col+1;     row24 = row+1; col24 = col+1;
%     %row31 = row-2; col31 = col;     row32 = row-1; col32 = col; 
%     %row34 = row+1; col34 = col;     row35 = row+2; col35 = col;
%     %row42 = row-1; col42 = col-1;   row43 = row; col43 = col-1;
%     %row44 = row+1; col44 = col-1;   row53 = row; col53 = col-2;
%     %row = [row;row13;row22;row23;row24;row31;row32;row34;row35;row42;row43;row44;row53];
%     %col = [col;col13;col22;col23;col24;col31;col32;col34;col35;col42;col43;col44;col53];
%      
%     BadptCol=[row]; BadptRow=[col]; 
%     [tempindex1] = find(BadptRow<1+max(coordinatesFEM(:,1)));  [tempindex2] = find(BadptCol<1+max(coordinatesFEM(:,2)));
%     [tempindex3] = find(BadptRow>min(coordinatesFEM(:,1)));  [tempindex4] = find(BadptCol>min(coordinatesFEM(:,2)));
%     tempindex1 = reshape(tempindex1,length(tempindex1),1);
%     tempindex2 = reshape(tempindex2,length(tempindex2),1);
%     tempindex3 = reshape(tempindex3,length(tempindex3),1);
%     tempindex4 = reshape(tempindex4,length(tempindex4),1);
%     tempindex = unique(intersect(tempindex4,intersect(tempindex3,intersect(tempindex1,tempindex2))));
%  
%     BadptRow = BadptRow(tempindex); BadptCol = BadptCol(tempindex); 
%     row=BadptRow; col=BadptCol;
%      
%      % --- Find which elements ---
%      BadptCoord = [];
%     for tempk = 1:length(row) % iterate for each clicking points
%         % Compute distance
%         DistCoordAll = (coordinatesFEM(:,2)-col(tempk)).^2 + (coordinatesFEM(:,1)-row(tempk)).^2; 
%         [DistCoordAllSorted, DistCoordAllSortedInd] = sort(sqrt(DistCoordAll));
%         ismember1 = []; ismember2 = []; ismember3 = []; 
%         for templ = 1:size(elementsFEM,1)
%             ismember1(templ) = ismember(DistCoordAllSortedInd(1),elementsFEM(templ,1:8));
%             ismember2(templ) = ismember(DistCoordAllSortedInd(2),elementsFEM(templ,1:8));
%             ismember3(templ) = ismember(DistCoordAllSortedInd(3),elementsFEM(templ,1:8));
%         end
%         [~,coll] = find(ismember1+ismember2+ismember3 == 3);
%         BadptCoord = [BadptCoord;elementsFEM(coll,:)'];
%     end
%     % Set unique BadptCoord
%     BadptCoord = unique(BadptCoord); BadptCoord = setdiff(BadptCoord,[0]);
%     
%     U(2*BadptCoord-1) = NaN; U(2*BadptCoord) = NaN; 
%     F(4*BadptCoord-3) = NaN; F(4*BadptCoord-2) = NaN; F(4*BadptCoord-1) = NaN; F(4*BadptCoord) = NaN;
% 
%     % ------ inpaint nans using gridfit ------
%     Coordxnodes = unique(coordinatesFEM(:,1)); Coordynodes = unique(coordinatesFEM(:,2));
%     nanindex = find(isnan(U(1:2:end))==1); notnanindex = setdiff([1:1:size(coordinatesFEM,1)],nanindex);
% %     [u1temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), U(2*notnanindex-1),Coordxnodes,Coordynodes,'regularizer','springs'); u1temp = u1temp';
% %     [v1temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), U(2*notnanindex),Coordxnodes,Coordynodes,'regularizer','springs'); v1temp = v1temp';
% %     [F11temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-3),Coordxnodes,Coordynodes,'regularizer','springs'); F11temp = F11temp';
% %     [F21temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-2),Coordxnodes,Coordynodes,'regularizer','springs'); F21temp = F21temp';
% %     [F12temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-1),Coordxnodes,Coordynodes,'regularizer','springs'); F12temp = F12temp';
% %     [F22temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-0),Coordxnodes,Coordynodes,'regularizer','springs'); F22temp = F22temp';
% 
%     [CoordxnodesGrid,CoordynodesGrid] = ndgrid(Coordxnodes,Coordynodes);
%     Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), U(2*notnanindex-1));
%     [u1temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %u1temp = u1temp';
%     Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), U(2*notnanindex));
%     [v1temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %v1temp = v1temp';
%     Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-3));
%     [F11temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %F11temp = F11temp';
%     Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-2));
%     [F21temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %F21temp = F21temp';
%     Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-1));
%     [F12temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %F12temp = F12temp';
%     Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-0));
%     [F22temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %F22temp = F22temp';
%      
%     
%     for tempi = 1:size(coordinatesFEM,1)
%         [row1,col1] = find(Coordxnodes==coordinatesFEM(tempi,1));
%         [row2,col2] = find(Coordynodes==coordinatesFEM(tempi,2));
%         U(2*tempi-1) = u1temp(row1,row2);
%         U(2*tempi)   = v1temp(row1,row2);
%         F(4*tempi-3) = F11temp(row1,row2); F(4*tempi-2) = F21temp(row1,row2);
%         F(4*tempi-1) = F12temp(row1,row2); F(4*tempi) = F22temp(row1,row2);
%     end
% 
%     close all; %Plotdisp_show(U,coordinatesFEM,elementsFEM);
%     Ftemp = F(2:4:end);
%     Sqx = zeros(4,size(elementsFEM,1)); Sqy = zeros(4,size(elementsFEM,1)); Sqc = zeros(4,size(elementsFEM,1));
%     for j = 1:size(elementsFEM,1)
%         Sqx(1:4,j) = coordinatesFEM(elementsFEM(j,1:4),1);
%         Sqy(1:4,j) = coordinatesFEM(elementsFEM(j,1:4),2);
%         Sqc(1:4,j) = Ftemp(elementsFEM(j,1:4));
%     end
%     if size(elementsFEM,1)>2e4
%         patch(Sqx,Sqy,Sqc,'facecolor','interp','edgecolor','none');
%     else
%         patch(Sqx,Sqy,Sqc,'facecolor','interp');
%     end
%     view(2); axis tight; axis equal; colorbar; xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
%     set(gcf,'color','w'); colormap jet;
%     title('Click bad F21 points and press -Enter- after done','fontweight','normal')
%  
%     prompt = 'Do you point out more F21-bad points? (0-yes; 1-no) Input: ';
%     ClearBadInitialPointsOrNot = input(prompt);
%     
% end
% 
% 
% 
% 
% %% ========= F12
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Have a look at integer search
% % --------------------------------------
% close all; %Plotdisp_show(U,coordinatesFEM,elementsFEM);
% Ftemp = F(3:4:end);
% Sqx = zeros(4,size(elementsFEM,1)); Sqy = zeros(4,size(elementsFEM,1)); Sqc = zeros(4,size(elementsFEM,1));
% for j = 1:size(elementsFEM,1)
%     Sqx(1:4,j) = coordinatesFEM(elementsFEM(j,1:4),1);
%     Sqy(1:4,j) = coordinatesFEM(elementsFEM(j,1:4),2);
%     Sqc(1:4,j) = Ftemp(elementsFEM(j,1:4));
% end
% if size(elementsFEM,1)>2e4
%     patch(Sqx,Sqy,Sqc,'facecolor','interp','edgecolor','none'); 
% else
%     patch(Sqx,Sqy,Sqc,'facecolor','interp'); 
% end
% view(2); axis tight; axis equal; colorbar; xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
% set(gcf,'color','w'); colormap jet;
% title('Click bad F12 points and press -Enter- after done','fontweight','normal')
% % figure; surf(v); colorbar; view(2)
% % title('Displacement v','fontweight','normal')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fprintf('Do you clear bad points by directly pointing F12-bad points? (0-yes; 1-no)  \n')
% prompt = 'Input here: ';
% ClearBadInitialPointsOrNot = input(prompt);
% 
% while ClearBadInitialPointsOrNot == 0
%      
%     [row1, col1] = ginput; row = floor(col1); col = floor(row1); 
%     %row13 = row; col13 = col+2;     row22 = row-1; col22 = col+1;
%     %row23 = row; col23 = col+1;     row24 = row+1; col24 = col+1;
%     %row31 = row-2; col31 = col;     row32 = row-1; col32 = col; 
%     %row34 = row+1; col34 = col;     row35 = row+2; col35 = col;
%     %row42 = row-1; col42 = col-1;   row43 = row; col43 = col-1;
%     %row44 = row+1; col44 = col-1;   row53 = row; col53 = col-2;
%     %row = [row;row13;row22;row23;row24;row31;row32;row34;row35;row42;row43;row44;row53];
%     %col = [col;col13;col22;col23;col24;col31;col32;col34;col35;col42;col43;col44;col53];
%      
%     BadptCol=[row]; BadptRow=[col]; 
%     [tempindex1] = find(BadptRow<1+max(coordinatesFEM(:,1)));  [tempindex2] = find(BadptCol<1+max(coordinatesFEM(:,2)));
%     [tempindex3] = find(BadptRow>min(coordinatesFEM(:,1)));  [tempindex4] = find(BadptCol>min(coordinatesFEM(:,2)));
%     tempindex1 = reshape(tempindex1,length(tempindex1),1);
%     tempindex2 = reshape(tempindex2,length(tempindex2),1);
%     tempindex3 = reshape(tempindex3,length(tempindex3),1);
%     tempindex4 = reshape(tempindex4,length(tempindex4),1);
%     tempindex = unique(intersect(tempindex4,intersect(tempindex3,intersect(tempindex1,tempindex2))));
%  
%     BadptRow = BadptRow(tempindex); BadptCol = BadptCol(tempindex); 
%     row=BadptRow; col=BadptCol;
%      
%      % --- Find which elements ---
%      BadptCoord = [];
%     for tempk = 1:length(row) % iterate for each clicking points
%         % Compute distance
%         DistCoordAll = (coordinatesFEM(:,2)-col(tempk)).^2 + (coordinatesFEM(:,1)-row(tempk)).^2; 
%         [DistCoordAllSorted, DistCoordAllSortedInd] = sort(sqrt(DistCoordAll));
%         ismember1 = []; ismember2 = []; ismember3 = []; 
%         for templ = 1:size(elementsFEM,1)
%             ismember1(templ) = ismember(DistCoordAllSortedInd(1),elementsFEM(templ,1:8));
%             ismember2(templ) = ismember(DistCoordAllSortedInd(2),elementsFEM(templ,1:8));
%             ismember3(templ) = ismember(DistCoordAllSortedInd(3),elementsFEM(templ,1:8));
%         end
%         [~,coll] = find(ismember1+ismember2+ismember3 == 3);
%         BadptCoord = [BadptCoord;elementsFEM(coll,:)'];
%     end
%     % Set unique BadptCoord
%     BadptCoord = unique(BadptCoord); BadptCoord = setdiff(BadptCoord,[0]);
%     
%     U(2*BadptCoord-1) = NaN; U(2*BadptCoord) = NaN; 
%     F(4*BadptCoord-3) = NaN; F(4*BadptCoord-2) = NaN; F(4*BadptCoord-1) = NaN; F(4*BadptCoord) = NaN;
% 
%     % ------ inpaint nans using gridfit ------
%     Coordxnodes = unique(coordinatesFEM(:,1)); Coordynodes = unique(coordinatesFEM(:,2));
%     nanindex = find(isnan(U(1:2:end))==1); notnanindex = setdiff([1:1:size(coordinatesFEM,1)],nanindex);
% %     [u1temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), U(2*notnanindex-1),Coordxnodes,Coordynodes,'regularizer','springs'); u1temp = u1temp';
% %     [v1temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), U(2*notnanindex),Coordxnodes,Coordynodes,'regularizer','springs'); v1temp = v1temp';
% %     [F11temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-3),Coordxnodes,Coordynodes,'regularizer','springs'); F11temp = F11temp';
% %     [F21temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-2),Coordxnodes,Coordynodes,'regularizer','springs'); F21temp = F21temp';
% %     [F12temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-1),Coordxnodes,Coordynodes,'regularizer','springs'); F12temp = F12temp';
% %     [F22temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-0),Coordxnodes,Coordynodes,'regularizer','springs'); F22temp = F22temp';
% 
%     [CoordxnodesGrid,CoordynodesGrid] = ndgrid(Coordxnodes,Coordynodes);
%     Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), U(2*notnanindex-1));
%     [u1temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %u1temp = u1temp';
%     Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), U(2*notnanindex));
%     [v1temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %v1temp = v1temp';
%     Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-3));
%     [F11temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %F11temp = F11temp';
%     Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-2));
%     [F21temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %F21temp = F21temp';
%     Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-1));
%     [F12temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %F12temp = F12temp';
%     Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-0));
%     [F22temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %F22temp = F22temp';
%      
%     
%     for tempi = 1:size(coordinatesFEM,1)
%         [row1,col1] = find(Coordxnodes==coordinatesFEM(tempi,1));
%         [row2,col2] = find(Coordynodes==coordinatesFEM(tempi,2));
%         U(2*tempi-1) = u1temp(row1,row2);
%         U(2*tempi)   = v1temp(row1,row2);
%         F(4*tempi-3) = F11temp(row1,row2); F(4*tempi-2) = F21temp(row1,row2);
%         F(4*tempi-1) = F12temp(row1,row2); F(4*tempi) = F22temp(row1,row2);
%     end
% 
%     close all; %Plotdisp_show(U,coordinatesFEM,elementsFEM);
%     Ftemp = F(3:4:end);
%     Sqx = zeros(4,size(elementsFEM,1)); Sqy = zeros(4,size(elementsFEM,1)); Sqc = zeros(4,size(elementsFEM,1));
%     for j = 1:size(elementsFEM,1)
%         Sqx(1:4,j) = coordinatesFEM(elementsFEM(j,1:4),1);
%         Sqy(1:4,j) = coordinatesFEM(elementsFEM(j,1:4),2);
%         Sqc(1:4,j) = Ftemp(elementsFEM(j,1:4));
%     end
%     if size(elementsFEM,1)>2e4
%         patch(Sqx,Sqy,Sqc,'facecolor','interp','edgecolor','none');
%     else
%         patch(Sqx,Sqy,Sqc,'facecolor','interp');
%     end
%     view(2); axis tight; axis equal; colorbar; xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
%     set(gcf,'color','w'); colormap jet;
%     title('Click bad F12 points and press -Enter- after done','fontweight','normal')
%  
%     prompt = 'Do you point out more F12-bad points? (0-yes; 1-no) Input: ';
%     ClearBadInitialPointsOrNot = input(prompt);
%     
% end
%   
% 
% 
% 
% %% ========= F22
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Have a look at integer search
% % --------------------------------------
% close all; %Plotdisp_show(U,coordinatesFEM,elementsFEM);
% Ftemp = F(4:4:end);
% Sqx = zeros(4,size(elementsFEM,1)); Sqy = zeros(4,size(elementsFEM,1)); Sqc = zeros(4,size(elementsFEM,1));
% for j = 1:size(elementsFEM,1)
%     Sqx(1:4,j) = coordinatesFEM(elementsFEM(j,1:4),1);
%     Sqy(1:4,j) = coordinatesFEM(elementsFEM(j,1:4),2);
%     Sqc(1:4,j) = Ftemp(elementsFEM(j,1:4));
% end
% if size(elementsFEM,1)>2e4
%     patch(Sqx,Sqy,Sqc,'facecolor','interp','edgecolor','none'); 
% else
%     patch(Sqx,Sqy,Sqc,'facecolor','interp'); 
% end
% view(2); axis tight; axis equal; colorbar; xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
% set(gcf,'color','w'); colormap jet;
% title('Click bad F22 points and press -Enter- after done','fontweight','normal')
% % figure; surf(v); colorbar; view(2)
% % title('Displacement v','fontweight','normal')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fprintf('Do you clear bad points by directly pointing F22-bad points? (0-yes; 1-no)  \n')
% prompt = 'Input here: ';
% ClearBadInitialPointsOrNot = input(prompt);
% 
% while ClearBadInitialPointsOrNot == 0
%      
%     [row1, col1] = ginput; row = floor(col1); col = floor(row1); 
%     %row13 = row; col13 = col+2;     row22 = row-1; col22 = col+1;
%     %row23 = row; col23 = col+1;     row24 = row+1; col24 = col+1;
%     %row31 = row-2; col31 = col;     row32 = row-1; col32 = col; 
%     %row34 = row+1; col34 = col;     row35 = row+2; col35 = col;
%     %row42 = row-1; col42 = col-1;   row43 = row; col43 = col-1;
%     %row44 = row+1; col44 = col-1;   row53 = row; col53 = col-2;
%     %row = [row;row13;row22;row23;row24;row31;row32;row34;row35;row42;row43;row44;row53];
%     %col = [col;col13;col22;col23;col24;col31;col32;col34;col35;col42;col43;col44;col53];
%      
%     BadptCol=[row]; BadptRow=[col]; 
%     [tempindex1] = find(BadptRow<1+max(coordinatesFEM(:,1)));  [tempindex2] = find(BadptCol<1+max(coordinatesFEM(:,2)));
%     [tempindex3] = find(BadptRow>min(coordinatesFEM(:,1)));  [tempindex4] = find(BadptCol>min(coordinatesFEM(:,2)));
%     tempindex1 = reshape(tempindex1,length(tempindex1),1);
%     tempindex2 = reshape(tempindex2,length(tempindex2),1);
%     tempindex3 = reshape(tempindex3,length(tempindex3),1);
%     tempindex4 = reshape(tempindex4,length(tempindex4),1);
%     tempindex = unique(intersect(tempindex4,intersect(tempindex3,intersect(tempindex1,tempindex2))));
%  
%     BadptRow = BadptRow(tempindex); BadptCol = BadptCol(tempindex); 
%     row=BadptRow; col=BadptCol;
%      
%      % --- Find which elements ---
%      BadptCoord = [];
%     for tempk = 1:length(row) % iterate for each clicking points
%         % Compute distance
%         DistCoordAll = (coordinatesFEM(:,2)-col(tempk)).^2 + (coordinatesFEM(:,1)-row(tempk)).^2; 
%         [DistCoordAllSorted, DistCoordAllSortedInd] = sort(sqrt(DistCoordAll));
%         ismember1 = []; ismember2 = []; ismember3 = []; 
%         for templ = 1:size(elementsFEM,1)
%             ismember1(templ) = ismember(DistCoordAllSortedInd(1),elementsFEM(templ,1:8));
%             ismember2(templ) = ismember(DistCoordAllSortedInd(2),elementsFEM(templ,1:8));
%             ismember3(templ) = ismember(DistCoordAllSortedInd(3),elementsFEM(templ,1:8));
%         end
%         [~,coll] = find(ismember1+ismember2+ismember3 == 3);
%         BadptCoord = [BadptCoord;elementsFEM(coll,:)'];
%     end
%     % Set unique BadptCoord
%     BadptCoord = unique(BadptCoord); BadptCoord = setdiff(BadptCoord,[0]);
%     
%     U(2*BadptCoord-1) = NaN; U(2*BadptCoord) = NaN; 
%     F(4*BadptCoord-3) = NaN; F(4*BadptCoord-2) = NaN; F(4*BadptCoord-1) = NaN; F(4*BadptCoord) = NaN;
% 
%     % ------ inpaint nans using gridfit ------
%     Coordxnodes = unique(coordinatesFEM(:,1)); Coordynodes = unique(coordinatesFEM(:,2));
%     nanindex = find(isnan(U(1:2:end))==1); notnanindex = setdiff([1:1:size(coordinatesFEM,1)],nanindex);
% %     [u1temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), U(2*notnanindex-1),Coordxnodes,Coordynodes,'regularizer','springs'); u1temp = u1temp';
% %     [v1temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), U(2*notnanindex),Coordxnodes,Coordynodes,'regularizer','springs'); v1temp = v1temp';
% %     [F11temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-3),Coordxnodes,Coordynodes,'regularizer','springs'); F11temp = F11temp';
% %     [F21temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-2),Coordxnodes,Coordynodes,'regularizer','springs'); F21temp = F21temp';
% %     [F12temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-1),Coordxnodes,Coordynodes,'regularizer','springs'); F12temp = F12temp';
% %     [F22temp] = gridfit(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-0),Coordxnodes,Coordynodes,'regularizer','springs'); F22temp = F22temp';
% 
%     [CoordxnodesGrid,CoordynodesGrid] = ndgrid(Coordxnodes,Coordynodes);
%     Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), U(2*notnanindex-1));
%     [u1temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); % u1temp = u1temp';
%     Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), U(2*notnanindex));
%     [v1temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %v1temp = v1temp';
%     Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-3));
%     [F11temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid);% F11temp = F11temp';
%     Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-2));
%     [F21temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %F21temp = F21temp';
%     Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-1));
%     [F12temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %F12temp = F12temp';
%     Ftemp = scatteredInterpolant(coordinatesFEM(notnanindex,1), coordinatesFEM(notnanindex,2), F(4*notnanindex-0));
%     [F22temp] = Ftemp(CoordxnodesGrid,CoordynodesGrid); %F22temp = F22temp';
%      
%     
%     for tempi = 1:size(coordinatesFEM,1)
%         [row1,col1] = find(Coordxnodes==coordinatesFEM(tempi,1));
%         [row2,col2] = find(Coordynodes==coordinatesFEM(tempi,2));
%         U(2*tempi-1) = u1temp(row1,row2);
%         U(2*tempi)   = v1temp(row1,row2);
%         F(4*tempi-3) = F11temp(row1,row2); F(4*tempi-2) = F21temp(row1,row2);
%         F(4*tempi-1) = F12temp(row1,row2); F(4*tempi) = F22temp(row1,row2);
%     end
% 
%     close all; %Plotdisp_show(U,coordinatesFEM,elementsFEM);
%     Ftemp = F(4:4:end);
%     Sqx = zeros(4,size(elementsFEM,1)); Sqy = zeros(4,size(elementsFEM,1)); Sqc = zeros(4,size(elementsFEM,1));
%     for j = 1:size(elementsFEM,1)
%         Sqx(1:4,j) = coordinatesFEM(elementsFEM(j,1:4),1);
%         Sqy(1:4,j) = coordinatesFEM(elementsFEM(j,1:4),2);
%         Sqc(1:4,j) = Ftemp(elementsFEM(j,1:4));
%     end
%     if size(elementsFEM,1)>2e4
%         patch(Sqx,Sqy,Sqc,'facecolor','interp','edgecolor','none');
%     else
%         patch(Sqx,Sqy,Sqc,'facecolor','interp');
%     end
%     view(2); axis tight; axis equal; colorbar; xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
%     set(gcf,'color','w'); colormap jet;
%     title('Click bad F22 points and press -Enter- after done','fontweight','normal')
%  
%     prompt = 'Do you point out more F22-bad points? (0-yes; 1-no) Input: ';
%     ClearBadInitialPointsOrNot = input(prompt);
%     
% end



  
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




