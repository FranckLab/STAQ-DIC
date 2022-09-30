% Generate a quadtree mesh considering sample's complex geometry
%
% ----------------------------------------------
% References
% [1] J Yang, K Bhattacharya. Fast adaptive mesh augmented Lagrangian Digital Image
% Correlation. Under review. 
% [2] S Funken, A Schmidt. Adaptive mesh refinement in 2D: an efficient
% implementation in MATLAB. Comp. Meth. Appl. Math. 20:459-479, 2020.
% ----------------------------------------------
% Author: Jin Yang 
% Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
% Last time updated: 02/2020.
% ==============================================


disp('--- Update the generated quadtree mesh ---')
%% ====== Remove finite elements where there is a hole ======
coordinatesFEMQuadtree = DICmesh.coordinatesFEM;
elementsFEMQuadtree = DICmesh.elementsFEM(:,1:4);
irregular = zeros(0,3);

while 1 % Generate a Quadtree mesh
    [~,mark4] = funMarkEdge(coordinatesFEMQuadtree,elementsFEMQuadtree,Df.ImgRefMask,DICmesh.elementMinSize*2); % Don't delete "*2"
    mark4 = find(mark4);
    [coordinatesFEMQuadtree,elementsFEMQuadtree,irregular] = QrefineR(coordinatesFEMQuadtree,elementsFEMQuadtree,irregular,mark4);
    if isempty(mark4)
        break
    end
end
  

%% %%%%% Plot refined mesh %%%%%
% figure; patch('Faces', elementsFEMQuadtree(:,1:4), 'Vertices', coordinatesFEMQuadtree, 'Facecolor','none','linewidth',1)
% axis equal; axis tight; set(gca,'fontsize',18); set(gcf,'color','w'); box on;
% xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
% title('Quadtree mesh','Interpreter','latex');
% a = gca; a.TickLabelInterpreter = 'latex';

% Update the quadtree mesh to deal with hanging nodes
for tempj=1:size(irregular,1)
    [Lia, Locb] = ismember( [irregular(tempj,1:2)], elementsFEMQuadtree(:,1:2), 'rows' );
    if Lia>0, elementsFEMQuadtree(Locb,8)=irregular(tempj,3); end
    [Lia, Locb] = ismember( [irregular(tempj,1:2)], elementsFEMQuadtree(:,2:3), 'rows' );
    if Lia>0, elementsFEMQuadtree(Locb,5)=irregular(tempj,3); end
    [Lia, Locb] = ismember( [irregular(tempj,1:2)], elementsFEMQuadtree(:,3:4), 'rows' );
    if Lia>0, elementsFEMQuadtree(Locb,6)=irregular(tempj,3); end
    [Lia, Locb] = ismember( [irregular(tempj,1:2)], elementsFEMQuadtree(:,[4,1]), 'rows' );
    if Lia>0, elementsFEMQuadtree(Locb,7)=irregular(tempj,3); end
    
    [Lia, Locb] = ismember( [irregular(tempj,1:2)], elementsFEMQuadtree(:,[2,1]), 'rows' );
    if Lia>0, elementsFEMQuadtree(Locb,8)=irregular(tempj,3); end
    [Lia, Locb] = ismember( [irregular(tempj,1:2)], elementsFEMQuadtree(:,[3,2]), 'rows' );
    if Lia>0, elementsFEMQuadtree(Locb,5)=irregular(tempj,3); end
    [Lia, Locb] = ismember( [irregular(tempj,1:2)], elementsFEMQuadtree(:,[4,3]), 'rows' );
    if Lia>0, elementsFEMQuadtree(Locb,6)=irregular(tempj,3); end
    [Lia, Locb] = ismember( [irregular(tempj,1:2)], elementsFEMQuadtree(:,[1,4]), 'rows' );
    if Lia>0, elementsFEMQuadtree(Locb,7)=irregular(tempj,3); end
end
 
% Remove elements within the center hole
[markInside4,markOutside4] = funMarkInside(coordinatesFEMQuadtree,elementsFEMQuadtree,Df.ImgRefMask);
elementsFEMQuadtree = elementsFEMQuadtree(markOutside4,:);


%% Fine nodes near the edges

% %%%%% Old codes: Erode ref image mask, and to find elements near holes' edges,
% nhood = [1 1 1; 1 1 1; 1 1 1];
% ImgRefMaskErode = DICpara.ImgRefMask;
% for tempi = 1: floor(0.5*max([20,mean(DICpara.winsize),mean(DICpara.winstepsize)]))-1
%     ImgRefMaskErode = imerode(ImgRefMaskErode, nhood);
% end
% % figure, imshow(ImgRefMaskErode');
% [markEleHoleEdge4,markEleFarOutside4] = funMarkInside(coordinatesFEMQuadtree,elementsFEMQuadtree,ImgRefMaskErode);
   
% %%%%% New codes: Find elements which are refined %%%%%%
elementsFEMQuadtreeSize = sqrt( ( coordinatesFEMQuadtree(elementsFEMQuadtree(:,1),1) - coordinatesFEMQuadtree(elementsFEMQuadtree(:,3),1) ).^2 + ...
        ( coordinatesFEMQuadtree(elementsFEMQuadtree(:,1),2) - coordinatesFEMQuadtree(elementsFEMQuadtree(:,3),2) ).^2 );
[markEleRefine4,~] = find(elementsFEMQuadtreeSize < 0.99*sqrt(2)*max([DICpara.winstepsize,0*DICpara.winsize]));
 
% %%%%% New codes: Find elements near the boudary %%%%%%
xMin = min(DICmesh.coordinatesFEM(:,1)); xMax = max(DICmesh.coordinatesFEM(:,1));
yMin = min(DICmesh.coordinatesFEM(:,2)); yMax = max(DICmesh.coordinatesFEM(:,2));
[row1,col1] = find(coordinatesFEMQuadtree(elementsFEMQuadtree(:,1),1) < xMin+1.01*DICpara.winstepsize);
[row2,col2] = find(coordinatesFEMQuadtree(elementsFEMQuadtree(:,3),1) > xMax-1.01*DICpara.winstepsize);
[row3,col3] = find(coordinatesFEMQuadtree(elementsFEMQuadtree(:,1),2) < yMin+1.01*DICpara.winstepsize);
[row4,col4] = find(coordinatesFEMQuadtree(elementsFEMQuadtree(:,3),2) > yMax-1.01*DICpara.winstepsize);
 
markEleHoleEdge4 =  union(row4,union(row3,union(row2,union(row1,markEleRefine4))));
markCoordHoleEdge = unique(elementsFEMQuadtree(markEleHoleEdge4,:));
try
    if markCoordHoleEdge(1)==0, markCoordHoleEdge = markCoordHoleEdge(2:end); end
catch
end

%%%%%% New codes: Find elements near marked elements %%%%%%
for tempi = 1 : 2 % 2+(round( 32 / mean(DICpara.winstepsize) )^2)
    
    markEleHoleEdgeNeigh4 = zeros(size(elementsFEMQuadtree,1),1);
    for eleInd = 1:size(elementsFEMQuadtree,1)
        markEleHoleEdgeNeigh4(eleInd) = length(intersect(elementsFEMQuadtree(eleInd,:),markCoordHoleEdge));
    end
    [markEleHoleEdgeNeigh4,~] = find(markEleHoleEdgeNeigh4>0);
    %%%%%%%%%
    markCoordHoleEdge = unique(elementsFEMQuadtree(markEleHoleEdgeNeigh4,:)) ;
    try
        if markCoordHoleEdge(1) == 0, markCoordHoleEdge = markCoordHoleEdge(2:end); end
    catch
    end
    
end


% %%%%% Store data structure %%%%%
DICmesh.markCoordHoleEdge = markCoordHoleEdge;
DICmesh.dirichlet = DICmesh.markCoordHoleEdge;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%% Plot %%%%%
figure; 
patch('Faces', elementsFEMQuadtree(:,1:4), 'Vertices', coordinatesFEMQuadtree, 'Facecolor','white','linewidth',1);
patch('Faces', elementsFEMQuadtree(markEleHoleEdge4,1:4), 'Vertices', coordinatesFEMQuadtree, 'Facecolor','yellow','linewidth',1);
hold on; patch('Faces', elementsFEMQuadtree(markEleHoleEdgeNeigh4,1:4), 'Vertices', coordinatesFEMQuadtree, 'Facecolor','yellow','linewidth',1);
axis equal; axis tight; set(gca,'fontsize',18); set(gcf,'color','w'); box on;
xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
title('Quadtree mesh','Interpreter','latex');
a = gca; a.TickLabelInterpreter = 'latex';

lgd = legend('Quadtree mesh elements','Elements near the edge','interpreter','latex','location','northeastoutside');
set(lgd,'fontsize',13);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialize variable U for the generated quadtree mesh
F_dispu = scatteredInterpolant( DICmesh.coordinatesFEM(:,1),DICmesh.coordinatesFEM(:,2),USubpb2(1:2:end) );
F_dispv = scatteredInterpolant( DICmesh.coordinatesFEM(:,1),DICmesh.coordinatesFEM(:,2),USubpb2(2:2:end) );
F_F11 = scatteredInterpolant( DICmesh.coordinatesFEM(:,1),DICmesh.coordinatesFEM(:,2),FSubpb2(1:4:end) );
F_F21 = scatteredInterpolant( DICmesh.coordinatesFEM(:,1),DICmesh.coordinatesFEM(:,2),FSubpb2(2:4:end) );
F_F12 = scatteredInterpolant( DICmesh.coordinatesFEM(:,1),DICmesh.coordinatesFEM(:,2),FSubpb2(3:4:end) );
F_F22 = scatteredInterpolant( DICmesh.coordinatesFEM(:,1),DICmesh.coordinatesFEM(:,2),FSubpb2(4:4:end) );

F_vdual_1 = scatteredInterpolant( DICmesh.coordinatesFEM(:,1),DICmesh.coordinatesFEM(:,2),vdual(1:2:end) );
F_vdual_2 = scatteredInterpolant( DICmesh.coordinatesFEM(:,1),DICmesh.coordinatesFEM(:,2),vdual(2:2:end) );
F_udual_11 = scatteredInterpolant( DICmesh.coordinatesFEM(:,1),DICmesh.coordinatesFEM(:,2),udual(1:4:end) );
F_udual_21 = scatteredInterpolant( DICmesh.coordinatesFEM(:,1),DICmesh.coordinatesFEM(:,2),udual(2:4:end) );
F_udual_12 = scatteredInterpolant( DICmesh.coordinatesFEM(:,1),DICmesh.coordinatesFEM(:,2),udual(3:4:end) );
F_udual_22 = scatteredInterpolant( DICmesh.coordinatesFEM(:,1),DICmesh.coordinatesFEM(:,2),udual(4:4:end) );
 
USubpb2 = 0*coordinatesFEMQuadtree(:);
temp = F_dispu(coordinatesFEMQuadtree(:,1),coordinatesFEMQuadtree(:,2)); USubpb2(1:2:end)=temp(:);
temp = F_dispv(coordinatesFEMQuadtree(:,1),coordinatesFEMQuadtree(:,2)); USubpb2(2:2:end)=temp(:);

FSubpb2 = repmat(0*coordinatesFEMQuadtree(:),2,1);
temp = F_F11(coordinatesFEMQuadtree(:,1),coordinatesFEMQuadtree(:,2)); FSubpb2(1:4:end)=temp(:);
temp = F_F21(coordinatesFEMQuadtree(:,1),coordinatesFEMQuadtree(:,2)); FSubpb2(2:4:end)=temp(:);
temp = F_F12(coordinatesFEMQuadtree(:,1),coordinatesFEMQuadtree(:,2)); FSubpb2(3:4:end)=temp(:);
temp = F_F22(coordinatesFEMQuadtree(:,1),coordinatesFEMQuadtree(:,2)); FSubpb2(4:4:end)=temp(:);

vdual = 0*coordinatesFEMQuadtree(:);
temp = F_vdual_1(coordinatesFEMQuadtree(:,1),coordinatesFEMQuadtree(:,2)); vdual(1:2:end)=temp(:);
temp = F_vdual_2(coordinatesFEMQuadtree(:,1),coordinatesFEMQuadtree(:,2)); vdual(2:2:end)=temp(:);

udual = repmat(0*coordinatesFEMQuadtree(:),2,1);
temp = F_udual_11(coordinatesFEMQuadtree(:,1),coordinatesFEMQuadtree(:,2)); udual(1:4:end)=temp(:);
temp = F_udual_21(coordinatesFEMQuadtree(:,1),coordinatesFEMQuadtree(:,2)); udual(2:4:end)=temp(:);
temp = F_udual_12(coordinatesFEMQuadtree(:,1),coordinatesFEMQuadtree(:,2)); udual(3:4:end)=temp(:);
temp = F_udual_22(coordinatesFEMQuadtree(:,1),coordinatesFEMQuadtree(:,2)); udual(4:4:end)=temp(:);




Plotdisp_show(USubpb2,coordinatesFEMQuadtree,elementsFEMQuadtree(:,1:4),DICpara,'EdgeColor');

DICmesh.coordinatesFEM = coordinatesFEMQuadtree;
DICmesh.elementsFEM = elementsFEMQuadtree;
DICmesh.coordinatesFEMWorld = [DICmesh.coordinatesFEM(:,1),size(DICpara.ImgRefMask,2)+1-DICmesh.coordinatesFEM(:,2)];
 

ResultFEMesh{1+floor(fNormalizedNewIndex/DICpara.ImgSeqIncUnit)} = ... % To save first mesh info
    struct( 'coordinatesFEM',DICmesh.coordinatesFEM,'elementsFEM',DICmesh.elementsFEM, ...
    'winsize',DICpara.winsize,'winstepsize',DICpara.winstepsize,'gridxyROIRange',DICpara.gridxyROIRange, ...
    'coordinatesFEMWorld',DICmesh.coordinatesFEMWorld,'elementMinSize',DICmesh.elementMinSize,'markCoordHoleEdge',DICmesh.markCoordHoleEdge);

 

% ====== Clear temporary variables ======
clear  irregular C R h mark4 markOutside4 Lia Locb F_dispu F_dispv temp


% ===== Remove bad points =====
[USubpb2,FSubpb2] = funRemoveOutliersQuadtree(DICmesh,DICpara,USubpb2,FSubpb2);
% [vdual,udual] = funRemoveOutliersQuadtree(DICmesh,DICpara,vdual,udual);

