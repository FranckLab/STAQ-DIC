function PlotaccelQuadtreeMasks(U,A,coordinatesFEM,elementsFEM,CurrentImg,CurrentImgMask,...
    DICpara,ImgSeqNum,Rnew,CircleFitPar)
%PLOTACCELQUADTREEMASKS: to plot DIC solved displacement components  
%   PlotaccelQuadtreeMasks(U,V,coordinatesFEMWorld,elementsFEM,CurrentImg,DICpara)
% ----------------------------------------------
%
%   INPUT: U                    Displacement vector: 
%                               U = [Ux_node1, Uy_node1, Ux_node2, Uy_node2, ... , Ux_nodeN, Uy_nodeN]';
%          A                    Accelerator vector:
%                               A = [Ax_node1, Ay_node1, Ax_node2, Ay_node2, ... , Ax_nodeN, Ay_nodeN]';
%          coordinatesFEM       FE mesh coordinates
%          elementsFEM          FE mesh elements
%          CurrentImg           Current deformed image
%          DICpara              DIC paramters
%
%   OUTPUT: Plots of x-displacement field and y-displacement field.
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
% Last date modified: 2020.12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%% TODO: %%%%%%%%%%%%

% ====== GE10 First Alex_cav_dic dataset ======
% load('D:\JinYang\MATLAB\2D_ALDIC_v4.2\img_alex_cav\R_data.mat');

% ====== GE10 150um_Top_2Mfps_1 =========
% load('D:\JinYang\MATLAB\2D_ALDIC_v4.2\img_cav_dic\img_cav_70121_DICpaper_150um_Top_2Mfps_1\R_data.mat');
% Rnew([2:6]) = [];
% CircleFitPar(2:6,:) = [];

% ====== GE10 100um_Top_2Mfps_1 =========
% load('D:\JinYang\MATLAB\2D_ALDIC_v4.2\img_cav_70121_DICpaper_100um_Top_2Mfps_1\R_data.mat');
% Rnew([2:6]) = [];
% CircleFitPar(2:6,:) = [];

% ===== GE10 19_53_04 =======
% load('D:\JinYang\MATLAB\2D_ALDIC_v4.2\img_19_53_04\R_data.mat');



% ===== GE4shot1 ========
% load('D:\JinYang\MATLAB\2D_ALDIC_v4.2\img_cav_dic\img_cav_92121_4percent_shot1\R_data.mat');
% Rnew([2:7]) = [];
% CircleFitPar(2:7,:) = [];

% ===== GE6shot3 12_58_29 ========
% load('D:\JinYang\MATLAB\2D_ALDIC_v4.2\img_cav_dic\img_cav_92121_6percent_shot3\R_data.mat');
% Rnew([2:7]) = [];
% CircleFitPar(2:7,:) = [];

% ===== GE10shot1 ========
% load('D:\JinYang\MATLAB\2D_ALDIC_v4.2\img_cav_dic\img_cav_92121_10percent_shot1\R_data.mat');
% Rnew([2:7]) = [];
% CircleFitPar(2:7,:) = [];

% ===== GE14shot1 ========
% load('D:\JinYang\MATLAB\2D_ALDIC_v4.2\img_cav_dic\img_cav_92121_14percent_shot1\R_data.mat');
% Rnew([2:7]) = [];
% CircleFitPar(2:7,:) = [];

% JY!!!!
% ImgSeqNum = 252-ImgSeqNum;
R = Rnew(ImgSeqNum);


%% Initialization
warning off; load('./plotFiles/colormap_RdYlBu.mat','cMap');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% convert pixel unit to the physical world unit %%%%%
try um2px = DICpara.um2px; 
catch um2px = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OrigDICImgTransparency = DICpara.OrigDICImgTransparency; % Original raw DIC image transparency
Image2PlotResults = DICpara.Image2PlotResults; % Choose image to plot over (first only, second and next images)

disp_u = U(1:2:end); disp_v = U(2:2:end);
accel_x = A(1:2:end); accel_y = A(2:2:end);

coordinatesFEMWorldDef = [coordinatesFEM(:,1)+Image2PlotResults*disp_u, coordinatesFEM(:,2)+Image2PlotResults*disp_v];

%%%%%%%%%%% JY!!!Mask START %%%%%%%%%%%%%%%
% if Image2PlotResults == 1
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
% end
%%%%%%%%%%% JY!!!Mask START %%%%%%%%%%%%%%%

 

bubble_y = (251-mean(CircleFitPar(end-20:end,1)))*DICpara.um2px;
bubble_x = mean(CircleFitPar(end-20:end,2))*DICpara.um2px;
% bubble_y = (251-28.725)*um2px;
% bubble_x = 194.33*um2px;

r = sqrt( (coordinatesFEMWorldDef(:,1)-bubble_x).^2 + (coordinatesFEMWorldDef(:,2)-bubble_y).^2  );
theta = atan2( (coordinatesFEMWorldDef(:,2)-bubble_y), coordinatesFEMWorldDef(:,1)-bubble_x);

% disp_rt = [cos(theta), -sin(theta); sin(theta), cos(theta)] * [disp_u, disp_v]';
% disp_r = disp_rt(1);
% disp_t = disp_rt(2);
disp_r = cos(theta).*disp_u + sin(theta).*disp_v; % JY!!! correction on 08/15/2021
disp_t = - sin(theta).*disp_u + cos(theta).*disp_v;
accel_r = cos(theta).*accel_x + sin(theta).*accel_y;
accel_t = - sin(theta).*accel_x + cos(theta).*accel_y;

[temp1,~] = find( r <  ( 0.0*DICpara.winsize + R + 0.1 )*DICpara.um2px );
disp_r( temp1 ) = nan;
disp_t( temp1 ) = nan;
accel_r( temp1 ) = nan;
accel_t( temp1 ) = nan;



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 1) dispx r ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1=figure; ax1=axes;  
try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
catch h1=surf(  flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis equal; axis tight; box on; set(gca,'fontSize',18); view(2);  set(gca,'ydir','normal');
hold on; ax2=axes; h2=show([],elementsFEM(:,1:4),coordinatesFEMWorldDef/um2px,accel_x,'NoEdgeColor');
set(gca,'fontSize',18); view(2); box on; axis equal;  axis tight;   
alpha(h2,OrigDICImgTransparency); colormap(cMap); caxis auto;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(cMap);  caxis([-10,40]); % colormap(jet);  
colormap(turbo); % caxis([0,70]); % caxis([-0.025,0.025]); 
% caxis([-1.3,-0.1]);
ax1.XTick = [100,200,300]; % Unit: px
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linkaxes([ax1,ax2]);  % Link axes together
ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
colormap(ax1,'gray'); % Give each one its own colormap
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
%%%%% convert pixel unit to the physical world unit %%%%%
xticklabels(ax1, num2cell(round(um2px*ax1.XTick*100)/100, length(ax1.XTick) )' );
yticklabels(ax1, num2cell(round(um2px*ax1.YTick*100)/100, length(ax1.YTick) )' );
cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';


 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 2) dispy t ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1=figure; ax1=axes;  
try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
catch h1=surf(  flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis equal; axis tight; box on; set(gca,'fontSize',18); view(2);  set(gca,'ydir','normal');
hold on; ax2=axes; h2=show([],elementsFEM(:,1:4),coordinatesFEMWorldDef/um2px,accel_y,'NoEdgeColor');
set(gca,'fontSize',18); view(2); box on; axis equal;  axis tight;   
alpha(h2,OrigDICImgTransparency); colormap(cMap); caxis auto;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(cMap ); % caxis([-2,2]);   % colormap(jet);  
% caxis([5,12]);
% caxis([0,38]);  
colormap(turbo); % caxis([-5,5]);
ax1.XTick = [100,200,300]; % Unit: px
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

linkaxes([ax1,ax2]);  % Link axes together
ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
colormap(ax1,'gray'); % Give each one its own colormap
set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
%%%%% convert pixel unit to the physical world unit %%%%%
xticklabels(ax1, num2cell(round(um2px*ax1.XTick*100)/100, length(ax1.XTick) )' );
yticklabels(ax1, num2cell(round(um2px*ax1.YTick*100)/100, length(ax1.YTick) )' );
cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';




