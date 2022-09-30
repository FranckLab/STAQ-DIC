function [strain_exx,strain_exy,strain_eyy,strain_principal_max,strain_principal_min, ...
    strain_maxshear,strain_vonMises] = PlotstrainQuadtreeMasks(U,F,coordinatesFEMWorld,elementsFEM,CurrentImg,CurrentImgMask,DICpara)
%PLOTSTRAINQUADTREE: to plot DIC solved strain fields on a quadtree mesh 
% and overlaid with the original DIC images 
%   [strain_exx,strain_exy,strain_eyy,strain_principal_max,strain_principal_min, ...
%    strain_maxshear,strain_vonMises] = PlotstrainQuadtree(U,F,coordinatesFEMWorld,elementsFEM,CurrentImg,DICpara)
% ----------------------------------------------
%
%   INPUT: F                    DIC solved deformation gradient tensor
%          coordinatesFE        FE mesh coordinates
%          elementsFEM          FE mesh elements
%
%   OUTPUT: strain_exx              strain xx-compoent
%           strain_exy              strain xy-compoent
%           strain_eyy              strain yy-compoent
%           strain_principal_max    max principal strain on the xy-plane
%           strain_principal_min    min principal strain on the xy-plane
%           strain_maxshear         max shear strain on the xy-plane
%           strain_vonMises         equivalent von Mises strain
%
%   Plots:       
%       1) strain sxx
%       2) strain sxy
%       3) strain syy
%       4) max principal strain on the xy-plane 
%       5) min principal strain on the xy-plane
%       6) max shear strain on the xy-plane
%       7) equivalent von Mises strain
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
warning off; load('./plotFiles/colormap_RdYlBu.mat','cMap');
% run('./img_Vito/Black_rainbow.m');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% convert pixel unit to the physical world unit %%%%%
try um2px = DICpara.um2px; 
catch um2px = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OrigDICImgTransparency = DICpara.OrigDICImgTransparency; % Original raw DIC image transparency
Image2PlotResults = DICpara.Image2PlotResults; % Choose image to plot over (first only, second and next images)

disp_u = U(1:2:end); disp_v = U(2:2:end);
coordinatesFEMWorldDef = [coordinatesFEMWorld(:,1)+Image2PlotResults*disp_u, coordinatesFEMWorld(:,2)+Image2PlotResults*disp_v];

%%%%%%%%%%% JY!!!Mask START %%%%%%%%%%%%%%%
if Image2PlotResults == 1
    if ~isempty(CurrentImgMask)
        for tempi = 1:size(coordinatesFEMWorldDef,1)
            try
            if CurrentImgMask( floor(coordinatesFEMWorldDef(tempi,1)/um2px), ...
                                (size(CurrentImgMask,2)+1-ceil(coordinatesFEMWorldDef(tempi,2)/um2px)) ) == 0 
                coordinatesFEMWorldDef(tempi,:) = [nan,nan];
            end
            catch
                coordinatesFEMWorldDef(tempi,:) = [nan,nan];
            end
        end
    else
        CurrentImgMask = imread(CurrentImg)';
        for tempi = 1:size(coordinatesFEMWorldDef,1)
            try
                if CurrentImgMask( floor(coordinatesFEMWorldDef(tempi,1)/um2px), ...
                        (size(CurrentImgMask,2)+1-ceil(coordinatesFEMWorldDef(tempi,2)/um2px)) ) < 0
                    coordinatesFEMWorldDef(tempi,:) = [nan,nan];
                end
            catch
                coordinatesFEMWorldDef(tempi,:) = [nan,nan];
            end
        end
        % [row1,~] = find(isnan(U(1:2:end))==1);
        % F(4*row1-3) = nan; F(4*row1-2) = nan; F(4*row1-1) = nan; F(4*row1) = nan; 
    end
end
%%%%%%%%%%% JY!!!Mask END %%%%%%%%%%%%%%%


%% Compute strain components

u_x = F(1:4:end); v_x = F(2:4:end);
u_y = F(3:4:end); v_y = F(4:4:end);

strain_exx = u_x; 
strain_exy = 0.5*(v_x+u_y);
strain_eyy = v_y;

strain_maxshear = sqrt((0.5*(strain_exx-strain_eyy)).^2 + strain_exy.^2);
% Principal strain
strain_principal_max = 0.5*(strain_exx+strain_eyy) + strain_maxshear;
strain_principal_min = 0.5*(strain_exx+strain_eyy) - strain_maxshear;
% equivalent von Mises strain
strain_vonMises = sqrt(strain_principal_max.^2 + strain_principal_min.^2 - ...
             strain_principal_max.*strain_principal_min + 3*strain_maxshear.^2);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 1) Strain exx ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1=figure; ax1=axes; 
try h1=imshow( flipud(imread(CurrentImg) ),'InitialMagnification','fit');
catch h1=surf( flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2=show([],elementsFEM(:,1:4),coordinatesFEMWorldDef/um2px,strain_exx,'NoEdgeColor');
set(gca,'fontSize',18); view(2); box on; axis equal;  axis tight;   
alpha(h2,OrigDICImgTransparency); colormap(turbo); caxis auto;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet); caxis([0,0.5]) % D Sample 
% colormap(jet);   caxis([-0.25,0.25]) % foam
% colormap(jet); caxis([-0.004,0]); % Sample 12
colormap(turbo); caxis([-0.5, 0.5])
% colormap(black_rainbow); caxis([-0.004,0.004]);
ax1.XTick = [100,200,300]; % Unit: px
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linkaxes([ax1,ax2]);  % Link axes together
ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
colormap(ax1,'gray'); % Give each one its own colormap
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
%%%%% convert pixel unit to the physical world unit %%%%%
xticklabels(ax1, num2cell(round(um2px*ax1.XTick*10)/10, length(ax1.XTick) )' );
yticklabels(ax1, num2cell(round(um2px*ax1.YTick*10)/10, length(ax1.YTick) )' );
cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';



 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 2) Strain exy ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1=figure; ax1=axes; 
try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
catch h1=surf( flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2=show([],elementsFEM(:,1:4),coordinatesFEMWorldDef/um2px,strain_exy,'NoEdgeColor');
set(gca,'fontSize',18); view(2); box on; axis equal;  axis tight;   
alpha(h2,OrigDICImgTransparency); colormap(turbo); caxis auto;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet); caxis([-0.08,0.08]) % D Sample 
% colormap(jet);  caxis([-0.25,0.25]) % foam
% colormap(jet); caxis([-0.008,0.008]); % Sample 12 
colormap(turbo); caxis([-0.8, 0.8])
%   colormap(black_rainbow); caxis([-0.0018,0.0018]);
ax1.XTick = [100,200,300]; % Unit: px
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linkaxes([ax1,ax2]);  % Link axes together
ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
colormap(ax1,'gray'); % Give each one its own colormap
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
%%%%% convert pixel unit to the physical world unit %%%%%
xticklabels(ax1, num2cell(round(um2px*ax1.XTick*10)/10, length(ax1.XTick) )' );
yticklabels(ax1, num2cell(round(um2px*ax1.YTick*10)/10, length(ax1.YTick) )' );
cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';



 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 3) Strain eyy ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1=figure; ax1=axes; 
try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
catch h1=surf( flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2=show([],elementsFEM(:,1:4),coordinatesFEMWorldDef/um2px,strain_eyy,'NoEdgeColor');
set(gca,'fontSize',18); view(2); box on; axis equal;  axis tight;   
alpha(h2,OrigDICImgTransparency); colormap(turbo); caxis auto;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet); caxis([-0.15,0]) % D Sample 
% colormap(jet);  caxis([-0.25,0.25])% foam
% colormap(jet); caxis([-0.002,0.017]); % Sample 12 
colormap(turbo);   caxis([-0.5,0.5])
%  colormap(black_rainbow); caxis([-0.0021,0.0021]);
ax1.XTick = [100,200,300]; % Unit: px
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linkaxes([ax1,ax2]);  % Link axes together
ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
colormap(ax1,'gray'); % Give each one its own colormap
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
%%%%% convert pixel unit to the physical world unit %%%%%
xticklabels(ax1, num2cell(round(um2px*ax1.XTick*10)/10, length(ax1.XTick) )' );
yticklabels(ax1, num2cell(round(um2px*ax1.YTick*10)/10, length(ax1.YTick) )' );
cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';

 

  
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % ====== 4) Strain e_principal_max ======
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fig1=figure; ax1=axes; 
% try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
% catch h1=surf(  flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
% end
% 
% axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
% hold on; ax2=axes; h2=show([],elementsFEM(:,1:4),coordinatesFEMWorldDef/um2px,strain_principal_max,'NoEdgeColor');
% set(gca,'fontSize',18); view(2); box on; axis equal;  axis tight;   
% alpha(h2,OrigDICImgTransparency); colormap(jet); caxis auto;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%% TODO: manually modify colormap and caxis %%%%%%
% % colormap(jet);  caxis auto; % D Sample 
% colormap(jet); caxis auto % foam
% % colormap(jet); caxis([0,0.02]); % Sample 12 
% % colormap(jet); caxis([-0.2,0.2])
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% linkaxes([ax1,ax2]);  % Link axes together
% ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
% colormap(ax1,'gray'); % Give each one its own colormap
% set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
% ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
% %%%%% convert pixel unit to the physical world unit %%%%%
% xticklabels(ax1, num2cell(round(um2px*ax1.XTick*10)/10, length(ax1.XTick) )' );
% yticklabels(ax1, num2cell(round(um2px*ax1.YTick*10)/10, length(ax1.YTick) )' );
% cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';
% 
% 
% 
%   
% % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % ====== 5) Strain e_principal_min ======
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fig1=figure; ax1=axes; 
% try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
% catch h1=surf(  flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
% end
% 
% axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
% hold on; ax2=axes; h2=show([],elementsFEM(:,1:4),coordinatesFEMWorldDef/um2px,strain_principal_min,'NoEdgeColor');
% set(gca,'fontSize',18); view(2); box on; axis equal;  axis tight;   
% alpha(h2,OrigDICImgTransparency); colormap(jet); caxis auto;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%% TODO: manually modify colormap and caxis %%%%%%
% % colormap(jet);  caxis auto; % D Sample 
% colormap(jet); caxis auto % foam
% % colormap(jet); caxis([-0.008,0]); % Sample 12 
% % colormap(jet); caxis([-0.2,0.2])
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% linkaxes([ax1,ax2]);  % Link axes together
% ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
% colormap(ax1,'gray'); % Give each one its own colormap
% set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
% ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
% %%%%% convert pixel unit to the physical world unit %%%%%
% xticklabels(ax1, num2cell(round(um2px*ax1.XTick*10)/10, length(ax1.XTick) )' );
% yticklabels(ax1, num2cell(round(um2px*ax1.YTick*10)/10, length(ax1.YTick) )' );
% cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';
% 
% 
% 
%  
% % %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % ====== 6) Strain e_max_shear ======
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fig1=figure; ax1=axes; 
% try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
% catch h1=surf(  flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
% end
% 
% axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
% hold on; ax2=axes; h2=show([],elementsFEM(:,1:4),coordinatesFEMWorldDef/um2px,strain_maxshear,'NoEdgeColor');
% set(gca,'fontSize',18); view(2); box on; axis equal;  axis tight;   
% alpha(h2,OrigDICImgTransparency); colormap(jet); caxis auto;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%% TODO: manually modify colormap and caxis %%%%%%
% % colormap(jet);  caxis auto; % D Sample 
% colormap(jet); caxis auto % foam
% % colormap(jet); caxis([0,0.011]); % Sample 12 
% % colormap(jet); caxis([0,0.2])
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% linkaxes([ax1,ax2]);  % Link axes together
% ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
% colormap(ax1,'gray'); % Give each one its own colormap
% set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
% ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
% %%%%% convert pixel unit to the physical world unit %%%%%
% xticklabels(ax1, num2cell(round(um2px*ax1.XTick*10)/10, length(ax1.XTick) )' );
% yticklabels(ax1, num2cell(round(um2px*ax1.YTick*10)/10, length(ax1.YTick) )' );
% cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';
% 
%  
%  
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ====== 7) von Mises equivalent strain ======
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fig1=figure; ax1=axes; 
% try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
% catch h1=surf(  flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
% end
% 
% axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2);  set(gca,'ydir','normal');
% hold on; ax2=axes; h2=show([],elementsFEM(:,1:4),coordinatesFEMWorldDef/um2px,strain_vonMises,'NoEdgeColor');
% set(gca,'fontSize',18); view(2); box on; axis equal;  axis tight;   
% alpha(h2,OrigDICImgTransparency); colormap(jet); caxis auto;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%% TODO: manually modify colormap and caxis %%%%%%
% % colormap(jet);  caxis auto; % D Sample 
% colormap(jet); caxis([0,0.5]); % foam
% % colormap(jet); caxis([0,0.025]); % Sample 12 
% % colormap(jet); caxis([0,0.18])
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% linkaxes([ax1,ax2]);  % Link axes together
% ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
% colormap(ax1,'gray'); % Give each one its own colormap
% set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
% ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
% %%%%% convert pixel unit to the physical world unit %%%%%
% xticklabels(ax1, num2cell(round(um2px*ax1.XTick*10)/10, length(ax1.XTick) )' );
% yticklabels(ax1, num2cell(round(um2px*ax1.YTick*10)/10, length(ax1.YTick) )' );
% cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';
% 




end
 
 
