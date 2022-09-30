


%% Change this line to: "try mex -O ba_interp2.cpp; catch; end"
% [Comment 2]: 
To use a slower bi-cubic spline interpolation instead of ba_interp2 (bi-cubic), 
% % cd("./Splines_interp/lib_matlab"); CompileLib; cd("../../");  % This line is to mex bi-cubic spline interpolations
% % addpath("./Splines_interp/lib_matlab"); % dbstop if error % % Old version codes.



%%
% %%%%% Plot elements near the hole edge
% markEleHoleEdge4 = zeros(size(DICmesh.elementsFEM,1),1);
% for eleInd = 1:size(DICmesh.elementsFEM,1)
%     markEleHoleEdge4(eleInd) = length(intersect(DICmesh.elementsFEM(eleInd,:),DICmesh.markCoordHoleEdge));
% end
% [markEleHoleEdge4,col5] = find(markEleHoleEdge4>0);
% 
% figure; patch('Faces', DICmesh.elementsFEM(markEleHoleEdge4,1:4), 'Vertices', DICmesh.coordinatesFEM, 'Facecolor','none','linewidth',1)
% axis equal; axis tight; set(gca,'fontsize',18); set(gcf,'color','w'); box on;
% xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
% title('Quadtree mesh','Interpreter','latex');
% a = gca; a.TickLabelInterpreter = 'latex';
 


%% Plot results for each ImgSeqNum
ImgSeqNum=27;
% for tempk=0:3, FSubpb2(4*DICmesh.markCoordHoleEdge-tempk) = nan; end

USubpb2 = ResultDisp{ImgSeqNum-1}.U;
FSubpb2 = ResultDefGrad{ImgSeqNum-1}.F;

coordinatesFEM = ResultFEMeshEachFrame{ImgSeqNum-1}.coordinatesFEM;
elementsFEM = ResultFEMeshEachFrame{ImgSeqNum-1}.elementsFEM;
try markCoordHoleEdge = ResultFEMeshEachFrame{ImgSeqNum-1}.markCoordHoleEdge; catch; end


% for tempk=0:3, FSubpb2(4*DICmesh.markCoordHoleEdge-tempk) = nan; end

USubpb2World = USubpb2; USubpb2World(2:2:end) = -USubpb2(2:2:end);
FSubpb2World = FSubpb2;
close all; Plotdisp_show(USubpb2World,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4));
Plotstrain_show(FSubpb2World,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4));

    

%% Plot quadtree
figure
Sqx = zeros(4,size(DICmesh.elementsFEM,1));
Sqy = zeros(4,size(DICmesh.elementsFEM,1));
Sqc = zeros(4,size(DICmesh.elementsFEM,1));
for j = 1:size(DICmesh.elementsFEM,1)
    Sqx(1:4,j) = DICmesh.coordinatesFEM(DICmesh.elementsFEM(j,1:4),1);
    Sqy(1:4,j) = DICmesh.coordinatesFEM(DICmesh.elementsFEM(j,1:4),2);
    Sqc(1:4,j) = ConvItPerEletemp(DICmesh.elementsFEM(j,1:4));
end
if size(DICmesh.elementsFEM,1)>2e4
    patch(Sqx,Sqy,Sqc,'facecolor','interp','edgecolor','none');
else
    patch(Sqx,Sqy,Sqc,'facecolor','interp');
end
view(2); axis tight; axis equal; colorbar; xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
set(gcf,'color','w'); colormap jet;

%%


BW = imread('text.png'); figure, imshow(BW);



