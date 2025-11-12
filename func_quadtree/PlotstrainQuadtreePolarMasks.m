function [strain_err,strain_ert,strain_ett,strain_principal_max,strain_principal_min, ...
    strain_maxshear,strain_vonMises] = PlotstrainQuadtreePolarMasks(U,F, ...
    coordinatesFEM,elementsFEM,CurrentImg,CurrentImgMask,DICpara,markCoordHoleStrain,ImgSeqNum, ...
    Rnew,CircleFitPar,t_exp,R_exp,Rdot_exp,Rdotdot_exp,Ginf,alpha_qKV,mu)
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

%%%%%%%%% TODO: %%%%%%%%%%%%

% ====== GE10 First Alex_cav_dic dataset ======
% load('D:\JinYang\MATLAB\2D_ALDIC_v4.2\img_alex_cav\R_data.mat');
% load('D:\JinYang\MATLAB\2D_ALDIC_v4.2\img_alex_cav\GE10_RofTdata_RtFit.mat')
% R = Rnew(ImgSeqNum); % unit: pixel

% ====== GE10 150um_Top_2Mfps_1 =========
% load('D:\JinYang\MATLAB\2D_ALDIC_v4.2\img_cav_dic\img_cav_70121_DICpaper_150um_Top_2Mfps_1\R_data.mat');
% Rnew([2:6]) = [];
% CircleFitPar(2:6,:) = [];
% R = Rnew(ImgSeqNum); % unit: pixel
% 
% load('D:\JinYang\MATLAB\2D_ALDIC_v4.2\img_cav_dic\img_cav_70121_DICpaper_150um_Top_2Mfps_1\GE10_RofTdata_RtFit.mat');
% R_exp([2:6]) = []; % unit: m
% Rdot_exp([2:6]) = []; % unit: m/s
% Rdotdot_exp([2:6]) = []; % unit: m/s^2
% t_exp([2:6]) = []; % unit: s

% ====== GE10 100um_Top_2Mfps_1 =========
% load('D:\JinYang\MATLAB\2D_ALDIC_v4.2\img_cav_70121_DICpaper_100um_Top_2Mfps_1\R_data.mat');
% Rnew([2:6]) = [];
% CircleFitPar(2:6,:) = [];
% R = Rnew(ImgSeqNum); % unit: pixel
% 
% load('D:\JinYang\MATLAB\2D_ALDIC_v4.2\img_cav_70121_DICpaper_100um_Top_2Mfps_1\GE10_RofTdata_RtFit.mat');
% R_exp([2:6]) = []; % unit: m
% Rdot_exp([2:6]) = []; % unit: m/s
% Rdotdot_exp([2:6]) = []; % unit: m/s^2
% t_exp([2:6]) = []; % unit: s


% ===== GE10 19_53_04 =======
% load('D:\JinYang\MATLAB\2D_ALDIC_v4.2\img_19_53_04\R_data.mat');
% load('D:\JinYang\MATLAB\2D_ALDIC_v4.2\img_19_53_04\GE10_RofTdata_RtFit.mat');
% R = Rnew(ImgSeqNum); % unit: pixel


% ===== GE4shot1 ========
% load('D:\JinYang\MATLAB\2D_ALDIC_v4.2\img_cav_dic\img_cav_92121_4percent_shot1\R_data.mat');
% Rnew([2:7]) = [];
% CircleFitPar(2:7,:) = [];
% R = Rnew(ImgSeqNum); % unit: pixel
% 
% load('D:\JinYang\MATLAB\2D_ALDIC_v4.2\img_cav_dic\img_cav_92121_4percent_shot1\GE_RofTdata_RtFit.mat');
% R_exp([2:7]) = []; % unit: m
% Rdot_exp([2:7]) = []; % unit: m/s
% Rdotdot_exp([2:7]) = []; % unit: m/s^2
% t_exp([2:7]) = []; % unit: s


% ===== GE6shot3 12_58_29 ========
% load('D:\JinYang\MATLAB\2D_ALDIC_v4.2\img_cav_dic\img_cav_92121_6percent_shot3\R_data.mat');
% Rnew([2:7]) = [];
% CircleFitPar(2:7,:) = [];
% R = Rnew(ImgSeqNum); % unit: pixel
% 
% load('D:\JinYang\MATLAB\2D_ALDIC_v4.2\img_cav_dic\img_cav_92121_6percent_shot3\GE_RofTdata_RtFit.mat');
% R_exp([2:7]) = []; % unit: m
% Rdot_exp([2:7]) = []; % unit: m/s
% Rdotdot_exp([2:7]) = []; % unit: m/s^2
% t_exp([2:7]) = []; % unit: s
% 
% %%%%% Surrounding material properties %%%%%%
% Ginf = 735.2; alpha_qKV= 9.05; mu = 0.13; %qKV model


% ===== GE10shot1 ========
% load('D:\JinYang\MATLAB\2D_ALDIC_v4.2\img_cav_dic\img_cav_92121_10percent_shot1\R_data.mat');
% Rnew([2:7]) = [];
% CircleFitPar(2:7,:) = [];
% R = Rnew(ImgSeqNum); % unit: pixel
% 
% load('D:\JinYang\MATLAB\2D_ALDIC_v4.2\img_cav_dic\img_cav_92121_10percent_shot1\GE_RofTdata_RtFit.mat');
% R_exp([2:7]) = []; % unit: m
% Rdot_exp([2:7]) = []; % unit: m/s
% Rdotdot_exp([2:7]) = []; % unit: m/s^2
% t_exp([2:7]) = []; % unit: s
% 
% %%%%% Surrounding material properties %%%%%%
% Ginf = 3.076e3; alpha_qKV = 5.01; mu = 0.153; %qKV model
 

% ===== GE14shot1 ========
% load('D:\JinYang\MATLAB\2D_ALDIC_v4.2\img_cav_dic\img_cav_92121_14percent_shot1\R_data.mat');
% Rnew([2:7]) = [];
% CircleFitPar(2:7,:) = [];
% R = Rnew(ImgSeqNum); % unit: pixel
%  
% 
% load('D:\JinYang\MATLAB\2D_ALDIC_v4.2\img_cav_dic\img_cav_92121_14percent_shot1\GE_RofTdata_RtFit.mat');
% R_exp([2:7]) = []; % unit: m
% Rdot_exp([2:7]) = []; % unit: m/s
% Rdotdot_exp([2:7]) = []; % unit: m/s^2
% t_exp([2:7]) = []; % unit: s
% 
% %%%%% Surrounding material properties %%%%%%
% Ginf = 6.713e3; alpha_qKV = 8.64; mu = 0.0982; %qKV model


%%%%% Surrounding material properties %%%%%%
% Ginf = 4.41e3; alpha_qKV = 4.582; mu = 0.0503; % qKV model
% Ginf = 45.03e3; mu = 0.089; alpha_qKV = 0; % NHKV model
rho = 1e3; % mass density
c_sound = 1506; % sound speed in gel

R = Rnew(ImgSeqNum);


%% Initialization
warning off; load('./plotFiles/colormap_RdYlBu.mat','cMap');
load('D:\MATLAB\2D_ALDIC_v4.2\plotFiles\cMap_viridis.mat');
load('D:\MATLAB\2D_ALDIC_v4.2\plotFiles\cMap_magma.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% convert pixel unit to the physical world unit %%%%%
um2px = DICpara.um2px; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OrigDICImgTransparency = DICpara.OrigDICImgTransparency; % Original raw DIC image transparency
Image2PlotResults = DICpara.Image2PlotResults; % Choose image to plot over (first only, second and next images)

disp_u = U(1:2:end); disp_v = U(2:2:end);
coordinatesFEMWorldDef = [coordinatesFEM(:,1)+Image2PlotResults*disp_u, coordinatesFEM(:,2)+Image2PlotResults*disp_v];

markCoordNotHoleStrain = setdiff( [1:1:size(coordinatesFEMWorldDef,1)], markCoordHoleStrain );


%%%%%%%%%%% JY!!!Mask START %%%%%%%%%%%%%%%
% if Image2PlotResults == 1
% for tempi = 1:size(coordinatesFEMWorldDef,1)
%     try
%     if CurrentImgMask( floor(coordinatesFEMWorldDef(tempi,1)/um2px), ...
%                         (size(CurrentImgMask,2)+1-floor(coordinatesFEMWorldDef(tempi,2)/um2px)) ) == 0 
%         coordinatesFEMWorldDef(tempi,:) = [nan,nan];
%     end
%     catch
%         coordinatesFEMWorldDef(tempi,:) = [nan,nan];
%     end
% end
% end
%%%%%%%%%%% JY!!!Mask START %%%%%%%%%%%%%%%





%% Convert to polar coordinate and compuate disp components

bubble_y = (251-mean(CircleFitPar(end-20:end,1)))*DICpara.um2px;
bubble_x = mean(CircleFitPar(end-20:end,2))*DICpara.um2px;

% bubble_y = (251-CircleFitPar(ImgSeqNum,1))*DICpara.um2px;
% bubble_x = CircleFitPar(ImgSeqNum,2)*DICpara.um2px;
% bubble_y = (251-28.725)*um2px;
% bubble_x = 194.33*um2px;

r = sqrt( (coordinatesFEMWorldDef(:,1)-bubble_x).^2 + (coordinatesFEMWorldDef(:,2)-bubble_y).^2  );
theta = atan2( (coordinatesFEMWorldDef(:,2)-bubble_y), coordinatesFEMWorldDef(:,1)-bubble_x);
cs = cos(theta).*sin(theta);
c2 = cos(theta).^2;
s2 = sin(theta).^2;


[temp1,~] = find(r <  ( 1*DICpara.winstepsize + 1 + R)*DICpara.um2px );
% [temp1,~] = find(r <  (  R )*DICpara.um2px );
F( 4*temp1-3 ) = nan;
F( 4*temp1-2 ) = nan;
F( 4*temp1-1 ) = nan;
F( 4*temp1  ) = nan;


%% Compute strain components

u_x = F(1:4:end); v_x = F(2:4:end);
u_y = F(3:4:end); v_y = F(4:4:end);

% Jacobian = (1+u_x).*(1+v_y) - v_x.*u_y;


strain_exx = u_x ; 
strain_exy = 0.5*(v_x+u_y);
strain_eyy = v_y;

strain_maxshear = sqrt((0.5*(strain_exx-strain_eyy)).^2 + strain_exy.^2);
% Principal strain
strain_principal_max = 0.5*(strain_exx+strain_eyy) + strain_maxshear;
strain_principal_min = 0.5*(strain_exx+strain_eyy) - strain_maxshear;
% equivalent von Mises strain
strain_vonMises = sqrt( strain_principal_max.^2 + strain_principal_min.^2 - ...
             strain_principal_max.*strain_principal_min ); 
% there was a mistake before:   previous  "+ 3*strain_maxshear.^2" term should not be included here.


%%%%% Strain polar components %%%%%%
% strain_err = u_x.*c2 - u_y.*cs - v_x.*cs + v_y.*s2;
% strain_ert = 0.5*(  u_x.*cs + u_y.*c2 - v_x.*s2 - v_y.*cs + ...
%                     u_x.*cs - u_y.*s2 + v_x.*c2 - v_y.*cs );
% strain_ett = u_x.*s2 + u_y.*cs + v_x.*cs + v_y.*c2;
strain_err =  u_x.*c2 + u_y.*cs + v_x.*cs + v_y.*s2 ;
strain_ert = 0.5*(  -u_x.*cs + u_y.*c2 - v_x.*s2 + v_y.*cs + ...
                    -u_x.*cs - u_y.*s2 + v_x.*c2 + v_y.*cs );
strain_ett =  u_x.*s2 - u_y.*cs - v_x.*cs + v_y.*c2; 

strain_logErr = log(1+strain_err);
strain_logEtt = log(1+strain_ett);
strain_logErt = log(1+strain_ert);

Jacobian = (1+strain_err).*(1+strain_ett).^2;

%%%%% Hydrostatic pressure %%%%%
F11 = strain_err + 1; 
lambda = sqrt(1./F11);
term1 = Ginf*(2.5 - 2./lambda + 0.5*lambda.^(-4) + alpha_qKV*(-177/20 + 0.75*lambda.^(-8) ...
                        - 0.4*lambda.^(-5) - 1.5*lambda.^(-4) + 6./lambda + 4*lambda ));
term2 = rho*( (R_exp(ImgSeqNum)^2*Rdotdot_exp(ImgSeqNum) + ...
        2*R_exp(ImgSeqNum)*Rdot_exp(ImgSeqNum)^2)./(r*1e-6) - ...
        0.5*R_exp(ImgSeqNum)^4*Rdot_exp(ImgSeqNum)^2./((r*1e-6).^4) );

I1 = F11.^2 + 2*lambda.^2;
term3 = Ginf/3 * (1 + alpha_qKV*(I1-3)) .* ( lambda.^(-4) + 2*lambda.^2 ); 
    
p_hydrostatic = term1 + term2 - term3;
p_hydrostatic_tilde = term1 + term2;

%%%%% Cauchy stress %%%%%
B11 = F11.^2; 
B22 = (strain_ett + 1).^2;
D11 = -2*R_exp(ImgSeqNum)^2*Rdot_exp(ImgSeqNum)./((r*1e-6).^3);
D22 = R_exp(ImgSeqNum)^2*Rdot_exp(ImgSeqNum)./((r*1e-6).^3);

stress_srr = Ginf*(1+alpha_qKV*(I1-3)).*B11 + 2*mu*D11 - p_hydrostatic_tilde;
stress_stt = Ginf*(1+alpha_qKV*(I1-3)).*B22 + 2*mu*D22 - p_hydrostatic_tilde;
stress_vonMises = sqrt( 0.5*(  2*(stress_srr-stress_stt).^2 + (0).^2 ) );

    
%%%%% Remove some bad points %%%%%%         
% strain_err( markCoordHoleStrain(:) ) = nan;      
% strain_ert( markCoordHoleStrain(:) ) = nan;  
% strain_ett( markCoordHoleStrain(:) ) = nan;  
% strain_maxshear( markCoordHoleStrain(:) ) = nan; 
% strain_principal_max( markCoordHoleStrain(:) ) = nan; 
% strain_principal_min( markCoordHoleStrain(:) ) = nan; 
% strain_vonMises( markCoordHoleStrain(:) ) = nan; 
% strain_logErr( markCoordHoleStrain(:) ) = nan;
         


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 1) Strain err ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fig1=figure; ax1=axes; 
% try h1=imshow( flipud( imread(CurrentImg) ),'InitialMagnification','fit');
% catch h1=surf( flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
% end
% 
% axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
% hold on; ax2=axes; h2=show([],elementsFEM(:,1:4),coordinatesFEMWorldDef/um2px, strain_err,'NoEdgeColor');
% set(gca,'fontSize',18); view(2); box on; axis equal;  axis tight;   
% alpha(h2,OrigDICImgTransparency); colormap(jet); caxis auto;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%% TODO: manually modify colormap and caxis %%%%%%
% % colormap(jet); caxis([0,0.5]) % D Sample 
% % colormap(jet); caxis([-0.1,0.02]) % foam
% % colormap(jet); caxis([-0.004,0]); % Sample 12
% colormap(cMap_viridis);  caxis([ -0.25 ,0 ]);  
% ax1.XTick = [100,200,300]; % Unit: px
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% linkaxes([ax1,ax2]);  % Link axes together
% ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
% colormap(ax1,'gray'); % Give each one its own colormap
% set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
% ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
% %%%%% convert pixel unit to the physical world unit %%%%%
% xticklabels(ax1, num2cell(round(um2px*ax1.XTick*100)/100, length(ax1.XTick) )' );
% yticklabels(ax1, num2cell(round(um2px*ax1.YTick*100)/100, length(ax1.YTick) )' );
% cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';



 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ====== 2) Strain ert ======
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1=figure; ax1=axes; 
try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
catch h1=surf( flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2=show([],elementsFEM(:,1:4),coordinatesFEMWorldDef/um2px,strain_logErt,'NoEdgeColor');
set(gca,'fontSize',18); view(2); box on; axis equal;  axis tight;   
alpha(h2,OrigDICImgTransparency); colormap(jet); caxis auto;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet); caxis([-0.08,0.08]) % D Sample 
% colormap(jet); caxis([-0.06,0.06]) % foam
% colormap(jet); caxis([-0.008,0.008]); % Sample 12 
colormap(jet);  caxis([-0.05,0.05]);
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

 
 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ====== 3) Strain logEtt ======
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1=figure; ax1=axes; 
try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
catch h1=surf( flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2=show([],elementsFEM(:,1:4),coordinatesFEMWorldDef/um2px,strain_logEtt,'NoEdgeColor');
set(gca,'fontSize',18); view(2); box on; axis equal;  axis tight;   
alpha(h2,OrigDICImgTransparency); colormap(jet); caxis auto;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet); caxis([-0.15,0]) % D Sample 
% colormap(jet); caxis([-0.05,0.2]) % foam
% colormap(jet); caxis([-0.002,0.017]); % Sample 12 
colormap(jet);  caxis([0, 0.3]);   % caxis([-0.1,0.1]);
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

 

  
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 4) Strain e_principal_max ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% % colormap(jet); caxis auto % foam
% % colormap(jet); caxis([0,0.02]); % Sample 12 
% colormap(jet); caxis([-0.2,0.2])
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% linkaxes([ax1,ax2]);  % Link axes together
% ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
% colormap(ax1,'gray'); % Give each one its own colormap
% set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
% ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
% %%%%% convert pixel unit to the physical world unit %%%%%
% xticklabels(ax1, num2cell(round(um2px*ax1.XTick*100)/100, length(ax1.XTick) )' );
% yticklabels(ax1, num2cell(round(um2px*ax1.YTick*100)/100, length(ax1.YTick) )' );
% cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';



  
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 5) Strain e_principal_min ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% % colormap(jet); caxis auto % foam
% % colormap(jet); caxis([-0.008,0]); % Sample 12 
% colormap(jet); caxis([-0.2,0.2])
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% linkaxes([ax1,ax2]);  % Link axes together
% ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
% colormap(ax1,'gray'); % Give each one its own colormap
% set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
% ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
% %%%%% convert pixel unit to the physical world unit %%%%%
% xticklabels(ax1, num2cell(round(um2px*ax1.XTick*100)/100, length(ax1.XTick) )' );
% yticklabels(ax1, num2cell(round(um2px*ax1.YTick*100)/100, length(ax1.YTick) )' );
% cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';



 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 6) Strain e_max_shear ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% % colormap(jet); caxis auto % foam
% % colormap(jet); caxis([0,0.011]); % Sample 12 
% colormap(jet); caxis([0,0.2])
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% linkaxes([ax1,ax2]);  % Link axes together
% ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
% colormap(ax1,'gray'); % Give each one its own colormap
% set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
% ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
% %%%%% convert pixel unit to the physical world unit %%%%%
% xticklabels(ax1, num2cell(round(um2px*ax1.XTick*100)/100, length(ax1.XTick) )' );
% yticklabels(ax1, num2cell(round(um2px*ax1.YTick*100)/100, length(ax1.YTick) )' );
% cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';

 
 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 7) von Mises equivalent strain ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fig1=figure; ax1=axes; 
% try h1=imshow( flipud(imread(CurrentImg)),'InitialMagnification','fit');
% catch h1=surf(  flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
% end
% 
% axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2);  set(gca,'ydir','normal');
% hold on; ax2=axes; h2=show([],elementsFEM(:,1:4), coordinatesFEMWorldDef/um2px,strain_vonMises,'NoEdgeColor');
% set(gca,'fontSize',18); view(2); box on; axis equal;  axis tight;   
% alpha(h2,OrigDICImgTransparency); colormap(jet); caxis auto;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%% TODO: manually modify colormap and caxis %%%%%%
% % colormap(jet);  caxis auto; % D Sample 
% % colormap(jet); caxis auto % foam
% % colormap(jet); caxis([0,0.025]); % Sample 12 
% colormap(cMap); % caxis([0,0.45])
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% linkaxes([ax1,ax2]);  % Link axes together
% ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
% colormap(ax1,'gray'); % Give each one its own colormap
% set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
% ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
% %%%%% convert pixel unit to the physical world unit %%%%%
% xticklabels(ax1, num2cell(round(um2px*ax1.XTick*100)/100, length(ax1.XTick) )' );
% yticklabels(ax1, num2cell(round(um2px*ax1.YTick*100)/100, length(ax1.YTick) )' );
% cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 8) Log strain Err ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1=figure; ax1=axes; 
try h1=imshow( flipud(imread(CurrentImg) ),'InitialMagnification','fit');
catch h1=surf( flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
end

axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
hold on; ax2=axes; h2=show([],elementsFEM(:,1:4),coordinatesFEMWorldDef/um2px, strain_logErr,'NoEdgeColor');
set(gca,'fontSize',18); view(2); box on; axis equal;  axis tight;   
alpha(h2,OrigDICImgTransparency); colormap(jet); caxis auto;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet); caxis([0,0.5]) % D Sample 
% colormap(jet); caxis([-0.1,0.02]) % foam
% colormap(jet); caxis([-0.004,0]); % Sample 12
colormap(jet);   caxis([-0.55, 0]);  
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

  
 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 9) Hydrostatic pressure ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fig1=figure; ax1=axes; 
% try h1=imshow( flipud( imread(CurrentImg) ),'InitialMagnification','fit');
% catch h1=surf( flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
% end
% 
% axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
% hold on; ax2=axes; h2=show([],elementsFEM(:,1:4),coordinatesFEMWorldDef/um2px, ...
%                          p_hydrostatic, 'NoEdgeColor');
%             % sign(p_hydrostatic).*log10( 1 + floor(abs(p_hydrostatic))), 'NoEdgeColor');
% set(gca,'fontSize',18); view(2); box on; axis equal;  axis tight;   
% alpha(h2,OrigDICImgTransparency); colormap(jet); caxis auto;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%% TODO: manually modify colormap and caxis %%%%%%
% % colormap(jet); caxis([0,0.5]) % D Sample 
% % colormap(jet); caxis([-0.1,0.02]) % foam
% % colormap(jet); caxis([-0.004,0]); % Sample 12
% colormap(jet);   caxis([-1e6,1e6]);
% % Black_rainbow; colormap(black_rainbow_plus); caxis([-4.5,4.5]);
% ax1.XTick = [100,200,300]; % Unit: px
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% linkaxes([ax1,ax2]);  % Link axes together
% ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
% colormap(ax1,'gray'); % Give each one its own colormap
% set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
% ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
% %%%%% convert pixel unit to the physical world unit %%%%%
% xticklabels(ax1, num2cell(round(um2px*ax1.XTick*100)/100, length(ax1.XTick) )' );
% yticklabels(ax1, num2cell(round(um2px*ax1.YTick*100)/100, length(ax1.YTick) )' );
% cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';





%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ====== 10) Cauchy conMises stress ======
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fig1=figure; ax1=axes; 
% try h1=imshow( flipud(imread(CurrentImg) ),'InitialMagnification','fit');
% catch h1=surf( flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
% end
% 
% axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
% hold on; ax2=axes; h2=show([],elementsFEM(:,1:4),coordinatesFEMWorldDef/um2px, ...
%     sign(stress_vonMises).*log10( 1 + floor(abs(stress_vonMises))),'NoEdgeColor');
% set(gca,'fontSize',18); view(2); box on; axis equal;  axis tight;   
% alpha(h2,OrigDICImgTransparency); colormap(jet); caxis auto;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%% TODO: manually modify colormap and caxis %%%%%%
% % colormap(jet); caxis([0,0.5]) % D Sample 
% % colormap(jet); caxis([-0.1,0.02]) % foam
% % colormap(jet); caxis([-0.004,0]); % Sample 12
% colormap(jet);   caxis([ 1.5  , 5.5 ]);
% % Black_rainbow; colormap(black_rainbow_plus); caxis([-4.5,4.5]);
% ax1.XTick = [100,200,300]; % Unit: px
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% linkaxes([ax1,ax2]);  % Link axes together
% ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
% colormap(ax1,'gray'); % Give each one its own colormap
% set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
% ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
% %%%%% convert pixel unit to the physical world unit %%%%%
% xticklabels(ax1, num2cell(round(um2px*ax1.XTick*100)/100, length(ax1.XTick) )' );
% yticklabels(ax1, num2cell(round(um2px*ax1.YTick*100)/100, length(ax1.YTick) )' );
% cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';
%%%%%%%%%%%%%%% comment end %%%%%%%%%%%%%
 


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 11) Cauchy stress stt ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fig1=figure; ax1=axes; 
% try h1=imshow( flipud(imread(CurrentImg) ),'InitialMagnification','fit');
% catch h1=surf( flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
% end
% 
% axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
% hold on; ax2=axes; h2=show([],elementsFEM(:,1:4),coordinatesFEMWorldDef/um2px, stress_stt,'NoEdgeColor');
% set(gca,'fontSize',18); view(2); box on; axis equal;  axis tight;   
% alpha(h2,OrigDICImgTransparency); colormap(jet); caxis auto;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%% TODO: manually modify colormap and caxis %%%%%%
% % colormap(jet); caxis([0,0.5]) % D Sample 
% % colormap(jet); caxis([-0.1,0.02]) % foam
% % colormap(jet); caxis([-0.004,0]); % Sample 12
% colormap(cMap); % caxis([-0.2,0.2]);
% ax1.XTick = [100,200,300]; % Unit: px
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% linkaxes([ax1,ax2]);  % Link axes together
% ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
% colormap(ax1,'gray'); % Give each one its own colormap
% set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
% ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
% %%%%% convert pixel unit to the physical world unit %%%%%
% xticklabels(ax1, num2cell(round(um2px*ax1.XTick*100)/100, length(ax1.XTick) )' );
% yticklabels(ax1, num2cell(round(um2px*ax1.YTick*100)/100, length(ax1.YTick) )' );
% cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 12) Jacobian ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fig1=figure; ax1=axes; 
% try h1=imshow( flipud(imread(CurrentImg) ),'InitialMagnification','fit');
% catch h1=surf( flipud( imread(CurrentImg) ),'EdgeColor','none','LineStyle','none');
% end
% 
% axis on; axis equal; axis tight; box on; set(gca,'fontSize',18); view(2); set(gca,'ydir','normal');
% hold on; ax2=axes; h2=show([],elementsFEM(:,1:4),coordinatesFEMWorldDef/um2px, Jacobian,'NoEdgeColor');
% set(gca,'fontSize',18); view(2); box on; axis equal;  axis tight;   
% alpha(h2,OrigDICImgTransparency); colormap(jet); caxis auto;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%% TODO: manually modify colormap and caxis %%%%%%
% % colormap(jet); caxis([0,0.5]) % D Sample 
% % colormap(jet); caxis([-0.1,0.02]) % foam
% % colormap(jet); caxis([-0.004,0]); % Sample 12
% colormap(jet); % caxis([0.85,1.3]); 
% ax1.XTick = [100,200,300]; % Unit: px
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% linkaxes([ax1,ax2]);  % Link axes together
% ax2.Visible = 'off'; ax2.XTick = []; ax2.YTick = []; % Hide the top axes
% colormap(ax1,'gray'); % Give each one its own colormap
% set([ax1,ax2],'Position',[.17 .11 .685 .815]);  
% ax1.Visible = 'on'; ax1.TickLabelInterpreter = 'latex'; 
% %%%%% convert pixel unit to the physical world unit %%%%%
% xticklabels(ax1, num2cell(round(um2px*ax1.XTick*100)/100, length(ax1.XTick) )' );
% yticklabels(ax1, num2cell(round(um2px*ax1.YTick*100)/100, length(ax1.YTick) )' );
% cb2 = colorbar('Position',[.17+0.685+0.012 .11 .03 .815]); cb2.TickLabelInterpreter = 'latex';


end
 
 
