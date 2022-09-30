% ---------------------------------------------
% Augmented Lagrangian Digital Image Correlation (ALDIC_Quadtree)
% using an adaptive quadtree mesh, which was automatically generated
% based on the DIC raw images
% 
% Author: Jin Yang, PhD @Caltech
% Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
% Date: 2015.04,06,07; 2016.03,04; 2020.11
% ---------------------------------------------

%% Section 1: Clear MATLAB environment & mex set up Spline interpolation  
close all; clear; clc; clearvars -global
fprintf('------------ Section 1 Start ------------ \n')
setenv('MW_MINGW64_LOC','C:\TDM-GCC-64');
try mex -O ba_interp2.cpp; catch; end  % mex set up ba_interp2.cpp script
% [Comment]: If this line reports error but it works before, 
% Change line 15 to: "try mex -O ba_interp2.cpp; catch; end"
addpath("./func",'./src','./plotFiles','./func_quadtree','./func_quadtree/refinement','./plotFiles/export_fig-d966721/'); 
addpath('./Images_Quadtree_demo/Images_Sample12'); % TODO: addpath("./YOUR IMAGE FOLDER"); 
addpath('./func/rbfinterp');
fprintf('------------ Section 1 Done ------------ \n \n')
 

%% Section 2: Load DIC parameters and set up DIC parameters 
fprintf('------------ Section 2 Start ------------ \n')
% ====== Read images ====== 
[file_name,Img,DICpara] = ReadImageQuadtree; % Load DIC raw images
DICpara.ImgSeqIncUnit = 1; % postprocess every two consecutive frames
DICpara.winsizeMin = 4; % the minimum element size in the adaptive quadtree mesh
    
% ====== Load mask files ======
[mask_file_name,ImgMask] = ReadImageMasks;  

% %%%%%% Uncomment lines below to change the DIC computing region (ROI) manually %%%%%%
% DICpara.gridxROIRange = [gridxROIRange1,gridxROIRange2]; DICpara.gridyROIRange = [Val1, Val2];
% E.g., gridxROIRange = [224,918]; gridyROIRange = [787,1162];

% ====== Normalize images: fNormalized = (f-f_avg)/(f_std) ======
[ImgNormalized,DICpara.gridxyROIRange] = funNormalizeImg(Img,DICpara.gridxyROIRange); 
 
% ====== Initialize variable storage ======
ResultDisp = cell(length(ImgNormalized)-1,1);    ResultDefGrad = cell(length(ImgNormalized)-1,1);
ResultStrain = cell(length(ImgNormalized)-1,1);  ResultStress = cell(length(ImgNormalized)-1,1);
ResultFEMeshEachFrame = cell(length(ImgNormalized)-1,1); % To store FE-mesh for each frame: needs future improvment to make it more efficient.
ResultFEMesh = cell(ceil((length(ImgNormalized)-1)/DICpara.ImgSeqIncUnit),1); % For incremental DIC mode
DICpara.SizeOfFFTSearchRegion = [8,8];
fprintf('------------ Section 2 Done ------------ \n \n')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To solve each frame in an image sequence
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ImgSeqNum =  2 : length(ImgNormalized) 

    close all; 
    DICpara.NewFFTSearch = 1; % Apply the new FFT search since the incremental disp is small between two consecutive frames 
    disp(['Current image frame #: ', num2str(ImgSeqNum),'/',num2str(length(ImgNormalized))]); % Report current frame #
    
    % ====== Compute image gradients ======
    fNormalizedMask = double( ImgMask{ImgSeqNum-1} ) ; % Load the mask file of previous frame
    fNormalized = ImgNormalized{ImgSeqNum-1} .* fNormalizedMask; % Load previous frame 
    Df = funImgGradient(fNormalized,fNormalized,fNormalizedMask); % Finite difference to compute image grayscale gradients;
    
    gNormalizedMask = double(ImgMask{ImgSeqNum}); % Load the mask file of current frame
    gNormalized = ImgNormalized{ ImgSeqNum } .* gNormalizedMask ; % Load current deformed image frame 
    
    DICpara.ImgRefMask = fNormalizedMask;
    
    figure, 
    subplot(2,2,1); imshow(fNormalized'); title('fNormalized'); colorbar;
    subplot(2,2,2); imshow(gNormalized'); title('gNormalized'); colorbar;
    subplot(2,2,3); imshow(fNormalizedMask'); title('f mask'); colorbar;
    subplot(2,2,4); imshow(gNormalizedMask'); title('g mask'); colorbar;
    
    
    %% Section 3: Compute an initial guess of the unknown displacement field
    fprintf('\n'); fprintf('------------ Section 3 Start ------------ \n')
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This section is to find or update an initial guess of the unknown displacements.
    % The key idea is to either to use a new FFT-based cross correlation peak fitting,
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ====== FFT-based cross correlation ======
    DICpara.InitFFTSearchMethod = 0;
    [DICpara,x0temp_f,y0temp_f,u_f,v_f,cc] = IntegerSearch(fNormalized,gNormalized,file_name,DICpara);
    
    %%%%% Interpolate to reference frame "f" coordinate system %%%%%
    xnodes = max([1+0.5*DICpara.winsize,DICpara.gridxyROIRange.gridx(1)])  ...
        : DICpara.winstepsize : min([size(fNormalized,1)-0.5*DICpara.winsize-1,DICpara.gridxyROIRange.gridx(2)]);
    ynodes = max([1+0.5*DICpara.winsize,DICpara.gridxyROIRange.gridy(1)])  ...
        : DICpara.winstepsize : min([size(fNormalized,2)-0.5*DICpara.winsize-1,DICpara.gridxyROIRange.gridy(2)]);
    
    [x0temp,y0temp] = ndgrid(xnodes,ynodes);
    u_f_NotNanInd = find(~isnan(u_f(:)));
    
    op1 = rbfcreate( [x0temp_f(u_f_NotNanInd),y0temp_f(u_f_NotNanInd)]',[u_f(u_f_NotNanInd)]','RBFFunction', 'thinplate');
    rbfcheck_maxdiff = rbfcheck(op1); % Check: rbf thin-plate interpolation  
    if rbfcheck_maxdiff > 1e-3, disp('Please check rbf interpolation! Pause here.'); pause; end
    u = rbfinterp([x0temp(:),y0temp(:)]', op1 );
    
    op2 = rbfcreate( [x0temp_f(u_f_NotNanInd),y0temp_f(u_f_NotNanInd)]',[v_f(u_f_NotNanInd)]','RBFFunction', 'thinplate');
    rbfcheck_maxdiff = rbfcheck(op2); % Check: rbf thin-plate interpolation  
    if rbfcheck_maxdiff > 1e-3, disp('Please check rbf interpolation! Pause here.'); pause; end
    v = rbfinterp([x0temp(:),y0temp(:)]', op2 );
    
    u = regularizeNd([x0temp(:),y0temp(:)],u(:),{xnodes',ynodes'},1e-3);
    v = regularizeNd([x0temp(:),y0temp(:)],v(:),{xnodes',ynodes'},1e-3);
        
    % ====== DIC uniform FE-mesh set up ======
    [DICmesh] = MeshSetUp(x0temp,y0temp,DICpara); % clear x0temp y0temp;
    % ====== Initial Value ======
    U0 = Init(u,v,cc.max,DICmesh.x0,DICmesh.y0,0); % [Temp code:] PlotuvInit;
    
    for tempi = 1:size(u,1)
        for tempj = 1:size(u,2)
            try
                if ~fNormalizedMask(x0temp(tempi,tempj),y0temp(tempi,tempj))
                    U0(2*(tempj+(tempi-1)*(size(u,2)))) = nan;
                    U0(2*(tempj+(tempi-1)*(size(u,2)))-1) = nan;
                end
            catch
            end
        end
    end
    
    
    % ====== Deal with incremental mode ======
    % %%%%% Old codes %%%%%
    fNormalizedNewIndex = ImgSeqNum-mod(ImgSeqNum-2,DICpara.ImgSeqIncUnit)-1;
    if DICpara.ImgSeqIncUnit == 1, fNormalizedNewIndex = fNormalizedNewIndex-1; end
    ResultFEMesh{1+floor(fNormalizedNewIndex/DICpara.ImgSeqIncUnit)} = ... % To save first mesh info
       struct( 'coordinatesFEM',DICmesh.coordinatesFEM,'elementsFEM',DICmesh.elementsFEM, ...
       'winsize',DICpara.winsize,'winstepsize',DICpara.winstepsize,'gridxyROIRange',DICpara.gridxyROIRange );
    
    % ====== Generate a quadtree mesh considering sample's complex geometry ======
    DICmesh.elementMinSize = DICpara.winsizeMin; % min element size in the refined quadtree mesh
    GenerateQuadtreeMesh; % Generate the quadtree mesh
    
    % ====== No need to update search region in the incremental mode ======
    % DICpara.SizeOfFFTSearchRegion = [ ceil( max( [max(3+abs(U0(1:2:end))), 3] ) ), ...
    %                                  ceil( max( [max(3+abs(U0(2:2:end))), 3] ) ) ];
    
    % ====== Stpre current mesh ======
    ResultFEMeshEachFrame{ImgSeqNum-1} = struct( 'coordinatesFEM',DICmesh.coordinatesFEM,'elementsFEM',DICmesh.elementsFEM,'markCoordHoleEdge',DICmesh.markCoordHoleEdge );
    fprintf('------------ Section 3 Done ------------ \n \n')

 
    %% Section 4: ALDIC Subproblem 1 -or- Local ICGN Subset DIC
    fprintf('------------ Section 4 Start ------------ \n')
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This section is to solve the first local step in ALDIC: Subproblem 1
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    % ====== ALStep 1 Subproblem1: Local Subset DIC ======
    mu=0; beta=0; tol=1e-2; ALSolveStep=1; ALSub1Time=zeros(6,1); ALSub2Time=zeros(6,1); 
    ConvItPerEle=zeros(size(DICmesh.coordinatesFEM,1),6); ALSub1BadPtNum=zeros(6,1);
    disp(['***** Start step',num2str(ALSolveStep),' Subproblem1 *****'])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ------ Start Local DIC IC-GN iteration ------
    [USubpb1,FSubpb1,HtempPar,ALSub1Timetemp,ConvItPerEletemp,LocalICGNBadPtNumtemp,markCoordHoleStrain] = ...
        LocalICGNQuadtree(U0,DICmesh.coordinatesFEM,Df,fNormalized,gNormalized,DICpara,'GaussNewton',tol);
    ALSub1Time(ALSolveStep) = ALSub1Timetemp; ConvItPerEle(:,ALSolveStep) = ConvItPerEletemp; ALSub1BadPtNum(ALSolveStep) = LocalICGNBadPtNumtemp; toc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ------  Manually find some bad points from Local Subset ICGN step ------
    % Comment these lines below if you don't have local bad points
    % %%%%% Comment START %%%%%%
    % USubpb1(2*DICmesh.markCoordHoleEdge-1:2*DICmesh.markCoordHoleEdge) = nan;
    % FSubpb1(4*DICmesh.markCoordHoleEdge-3:4*DICmesh.markCoordHoleEdge) = nan;
    % [USubpb1,FSubpb1] = funRemoveOutliersQuadtree(DICmesh,DICpara,USubpb1,FSubpb1);
    % disp('--- Remove bad points done ---')
    % %%%%% Comment END %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    % %%%%% Convergence iteration # %%%%%
    % Convtemp = [ConvItPerEletemp, ConvItPerEletemp]'; Convtemp = Convtemp(:);
    % Plotdisp_show(Convtemp,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4),DICpara,'EdgeColor');
    
    % ====== Thin-plate interpolate bad points =====
    coordinatesFEM = DICmesh.coordinatesFEM;
    U = USubpb1; F = FSubpb1;
    nanindexU = find(isnan(U(1:2:end))==1); notnanindex = setdiff([1:1:size(coordinatesFEM,1)],nanindexU);
    nanindexF = find(isnan(F(1:4:end))==1); notnanindexF = setdiff([1:1:size(coordinatesFEM,1)],nanindexF);
    
    %%%%% Ux %%%%%%
    op1 = rbfcreate( [coordinatesFEM(notnanindex,1:2)]',[U(2*notnanindex-1)]','RBFFunction', 'thinplate');
    rbfcheck_maxdiff = rbfcheck(op1); % check rbf interpolation
    if rbfcheck_maxdiff > 1e-3, disp('Please check rbf interpolation! Pause here.'); pause; end
    fi1 = rbfinterp([coordinatesFEM(:,1:2)]', op1 );
    
    %%%%% Uy %%%%%%
    op2 = rbfcreate( [coordinatesFEM(notnanindex,1:2)]',[U(2*notnanindex)]','RBFFunction', 'thinplate');
    rbfcheck_maxdiff = rbfcheck(op2); % check rbf interpolation
    if rbfcheck_maxdiff > 1e-3, disp('Please check rbf interpolation! Pause here.'); pause; end
    fi2 = rbfinterp([coordinatesFEM(:,1:2)]', op2 );
    
    %%%%% Assemble [Ux, Uy] %%%%%
    U_rbf_thinplate = [fi1(:),fi2(:)]';  U_rbf_thinplate = U_rbf_thinplate(:);
    
    %%%%% F11 %%%%%
    op =rbfcreate([coordinatesFEM(notnanindex,1:2)]', ...
        F(4*notnanindex-3)','RBFFunction', 'thinplate');
    rbfcheck_maxdiff = rbfcheck(op);
    if rbfcheck_maxdiff > 1e-3, disp('Please check rbf interpolation! Pause here.'); pause; end
    fi11 = rbfinterp([coordinatesFEM(:,1:2)]', op );
    
    %%%%% F21 %%%%%
    op =rbfcreate([coordinatesFEM(notnanindex,1:2)]', ...
        F(4*notnanindex-2)','RBFFunction', 'thinplate');
    rbfcheck_maxdiff = rbfcheck(op);
    if rbfcheck_maxdiff > 1e-3, disp('Please check rbf interpolation! Pause here.'); pause; end
    fi21 = rbfinterp([coordinatesFEM(:,1:2)]', op );
    
    %%%%% F12 %%%%%
    op =rbfcreate([coordinatesFEM(notnanindex,1:2)]', ...
        F(4*notnanindex-1)','RBFFunction', 'thinplate');
    rbfcheck_maxdiff = rbfcheck(op);
    if rbfcheck_maxdiff > 1e-3, disp('Please check rbf interpolation! Pause here.'); pause; end
    fi12 = rbfinterp([coordinatesFEM(:,1:2)]', op );
    
    %%%%% F22 %%%%%
    op =rbfcreate([coordinatesFEM(notnanindex,1:2)]', ...
        F(4*notnanindex-0)','RBFFunction', 'thinplate');
    rbfcheck_maxdiff = rbfcheck(op);
    if rbfcheck_maxdiff > 1e-3, disp('Please check rbf interpolation! Pause here.'); pause; end
    fi22 = rbfinterp([coordinatesFEM(:,1:2)]', op );
    
    %%%%% Assemble [F11,F21,F12,F22] %%%%%
    F_rbf_thinplate = [fi11(:),fi21(:),fi12(:),fi22(:)]';  F_rbf_thinplate = F_rbf_thinplate(:);
     
    % ------ Plot ------
    USubpb1 = U_rbf_thinplate; FSubpb1 = F_rbf_thinplate;
    USubpb1World = USubpb1; USubpb1World(2:2:end) = -USubpb1(2:2:end); FSubpb1World = FSubpb1;
    
    Plotdisp_show(USubpb1World,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4),DICpara,'EdgeColor');
    Plotstrain_show(FSubpb1World,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4),DICpara,'EdgeColor');
    save(['Subpb1_step',num2str(ALSolveStep)],'USubpb1','FSubpb1');
    fprintf('------------ Section 4 Done ------------ \n \n')
    
   
    %% Section 5: Subproblem 2 -- solve the global compatible displacement field
    fprintf('------------ Section 5 Start ------------ \n'); tic;
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This section is to solve the global step in ALDIC Subproblem 2
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    % ======= ALStep 1 Subproblem 2: Global constraint =======
    % ------ Smooth displacements for a better F ------
    DICpara.DispFilterSize=0; DICpara.DispFilterStd=0; DICpara.StrainFilterSize=0; DICpara.StrainFilterStd=0; LevelNo=1;
    DICpara.DispSmoothness = 0; DICpara.StrainSmoothness = 1e-4;
    if DICpara.DispSmoothness>1e-6, USubpb1 = funSmoothDispQuadtree(USubpb1,DICmesh,DICpara); end
    if DICpara.StrainSmoothness>1e-6, FSubpb1 = funSmoothStrainQuadtree(FSubpb1,DICmesh,DICpara); end
    
	% ====== Define penalty parameter ======
    mu = 1e-3; udual = 0*FSubpb1; vdual = 0*USubpb1; 
    betaList = [1e-3,1e-2,1e-1]*mean(DICpara.winstepsize).^2.*mu; % Tune beta in the betaList 
    Err1 = zeros(length(betaList),1); Err2 = Err1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(['***** Start step',num2str(ALSolveStep),' Subproblem2 *****']);
    DICpara.GaussPtOrder = 2; alpha = 0;  % No regularization added
    % ====== Solver using finite element method ======
    if ImgSeqNum == 2
        for tempk = 1:length(betaList)
            beta = betaList(tempk); display(['Try #',num2str(tempk),' beta = ',num2str(beta)]);
            GaussPtOrder=3; alpha=0; [USubpb2] = Subpb2Quadtree(DICmesh,DICpara.GaussPtOrder,beta,mu,USubpb1,FSubpb1,udual,vdual,alpha,mean(DICpara.winstepsize),0);
            [FSubpb2,~,~] = funGlobalNodalStrainQuadtree(DICmesh,USubpb2,DICpara.GaussPtOrder,0);

            Err1(tempk) = norm(USubpb1-USubpb2,2);
            Err2(tempk) = norm(FSubpb1-FSubpb2,2);
        end
        
        Err1Norm = (Err1-mean(Err1))/std(Err1); % figure, plot(Err1Norm);
        Err2Norm = (Err2-mean(Err2))/std(Err2); % figure, plot(Err2Norm);
        ErrSum = Err1Norm+Err2Norm; % figure, plot(ErrSum); title('Tune the best \beta in the subproblem 2'); 
        [~,indexOfbeta] = min(ErrSum); 
     
        try % Tune the best beta by a quadratic polynomial 0fitting
            [fitobj] = fit(log10(betaList(indexOfbeta-1:1:indexOfbeta+1))',ErrSum(indexOfbeta-1:1:indexOfbeta+1),'poly2');
            p = coeffvalues(fitobj); beta = 10^(-p(2)/2/p(1));
        catch, beta = betaList(indexOfbeta);
        end
        display(['Best beta = ',num2str(beta)]);
    else 
        try beta = DICpara.beta;
        catch, beta = 1e-3*mean(DICpara.winstepsize).^2.*mu;
        end
    end
      
    % Using the optimal beta to solve the ALDIC Subproblem 2 again
    if abs(beta-betaList(end))>abs(eps)
        [USubpb2] = Subpb2Quadtree(DICmesh,DICpara.GaussPtOrder,beta,mu,USubpb1,FSubpb1,udual,vdual,alpha,mean(DICpara.winstepsize),0);
        [FSubpb2,~,~] = funGlobalNodalStrainQuadtree(DICmesh,USubpb2,DICpara.GaussPtOrder,0);
        ALSub2Time(ALSolveStep) = toc; toc
    end
    
    % ------- Smooth strain field --------
    if DICpara.DispSmoothness>1e-6, USubpb2 = funSmoothDispQuadtree(USubpb2,DICmesh,DICpara); end
    % ------- Don't smooth strain fields near the boundary --------
    for tempk=0:3, FSubpb2(4*DICmesh.markCoordHoleEdge-tempk) = FSubpb1(4*DICmesh.markCoordHoleEdge-tempk); end
    if DICpara.StrainSmoothness>1e-6, FSubpb2 = funSmoothStrainQuadtree(0.1*FSubpb2+0.9*FSubpb1,DICmesh,DICpara); end
    for tempk=0:1, USubpb2(2*markCoordHoleStrain-tempk) = USubpb1(2*markCoordHoleStrain-tempk); end
    for tempk=0:3, FSubpb2(4*markCoordHoleStrain-tempk) = FSubpb1(4*markCoordHoleStrain-tempk); end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    % ------- Save data ------
    save(['Subpb2_step',num2str(ALSolveStep)],'USubpb2','FSubpb2');

    % ------ Plot ------
    USubpb2World = USubpb2; USubpb2World(2:2:end) = -USubpb2(2:2:end); FSubpb2World = FSubpb2;  
    % close all; Plotdisp_show(USubpb2World,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4),DICpara,'EdgeColor');
    % Plotstrain_show(FSubpb2World,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4),DICpara,'EdgeColor');

    % ======= Update dual variables =======
    udual = FSubpb2 - FSubpb1; vdual = USubpb2 - USubpb1;
	save(['uvdual_step',num2str(ALSolveStep)],'udual','vdual');
    fprintf('------------ Section 5 Done ------------ \n \n')
 

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Section 6: ADMM iterations
    fprintf('------------ Section 6 Start ------------ \n')
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This section is the ADMM iteration, where both Subproblems 1 & 2 are solved iteratively.
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % ==================== ADMM AL Loop ==========================
    ALSolveStep = 1; tol2 = 1e-2; UpdateY = 1e4;  
    HPar = cell(21,1); for tempj = 1:21, HPar{tempj} = HtempPar(:,tempj); end

    while (ALSolveStep < 3)
        ALSolveStep = ALSolveStep + 1;  % Update using the last step
        
        %%%%%%%% These lines can be used to further update each DIC subset window size %%%%%%%
        % Ftemp1 = FSubpb2(1:2:end); Ftemp2 = FSubpb2(2:2:end);
        % [DFtemp1,~,~] = funGlobalNodalStrainQuadtree(DICmesh,Ftemp1,DICpara.GaussPtOrder,0);
        % [DFtemp2,~,~] = funGlobalNodalStrainQuadtree(DICmesh,Ftemp2,DICpara.GaussPtOrder,0);
        %
        % winsize_x_ub1 = abs(2*FSubpb2(1:4:end)./DFtemp1(1:4:end)); winsize_x_ub2 = abs(2*FSubpb2(3:4:end)./DFtemp1(3:4:end));
        % winsize_y_ub1 = abs(2*FSubpb2(1:4:end)./DFtemp1(2:4:end)); winsize_y_ub2 = abs(2*FSubpb2(3:4:end)./DFtemp1(4:4:end));
        % winsize_x_ub3 = abs(2*FSubpb2(2:4:end)./DFtemp2(1:4:end)); winsize_x_ub4 = abs(2*FSubpb2(4:4:end)./DFtemp2(3:4:end));
        % winsize_y_ub3 = abs(2*FSubpb2(2:4:end)./DFtemp2(2:4:end)); winsize_y_ub4 = abs(2*FSubpb2(4:4:end)./DFtemp2(4:4:end));
        %
        % winsize_x_ub = round(min([winsize_x_ub1,winsize_x_ub2,winsize_x_ub3,winsize_x_ub4,DICpara.winsize*ones(length(winsize_x_ub1),1)],[],2));
        % winsize_x_List = max([winsize_x_ub, 10*ones(length(winsize_x_ub1),1)],[],2);
        % winsize_y_ub = round(min([winsize_y_ub1,winsize_y_ub2,winsize_y_ub3,winsize_y_ub4,DICpara.winsize*ones(length(winsize_y_ub1),1)],[],2));
        % winsize_y_List = max([winsize_y_ub, 10*ones(length(winsize_y_ub1),1)],[],2);
        % winsize_List = 2*ceil([winsize_x_List,winsize_y_List]/2);
        winsize_List = DICpara.winsize*ones(size(DICmesh.coordinatesFEM,1),2);
        DICpara.winsize_List = winsize_List;
        

        %%%%%%%%%%%%%%%%%%%%%%% Subproblem 1 %%%%%%%%%%%%%%%%%%%%%%%%%
        disp(['***** Start step',num2str(ALSolveStep),' Subproblem1 *****']);
        tic; [USubpb1,~,ALSub1Timetemp,ConvItPerEletemp,LocalICGNBadPtNumtemp] = Subpb1Quadtree(...
                                            USubpb2,FSubpb2,udual,vdual,DICmesh.coordinatesFEM,...
                                            Df,fNormalized,gNormalized,mu,beta,HPar,ALSolveStep,DICpara,'GaussNewton',tol);
        FSubpb1 = FSubpb2; toc 
        % for tempk=0:1, USubpb1(2*markCoordHoleStrain-tempk) = USubpb2(2*markCoordHoleStrain-tempk); end
        ALSub1Time(ALSolveStep) = ALSub1Timetemp; ConvItPerEle(:,ALSolveStep) = ConvItPerEletemp; ALSub1BadPtNum(ALSolveStep) = LocalICGNBadPtNumtemp;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ------  Manually find some bad points from Local Subset ICGN step ------
        % disp('--- Start to manually remove bad points --- \n')
        % disp('    Comment codes here if you do not have bad local points \n')
        % %%%%% Comment START %%%%%
        %  [USubpb1,FSubpb1] = funRemoveOutliersQuadtree(DICmesh,DICpara,USubpb1,FSubpb1);
        %  disp('--- Remove bad points done ---')
        % %%%%% Comment END %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        save(['Subpb1_step',num2str(ALSolveStep)],'USubpb1','FSubpb1');
         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ============== Subproblem 2 ==============
        disp(['***** Start step',num2str(ALSolveStep),' Subproblem2 *****'])
        tic; [USubpb2] = Subpb2Quadtree(DICmesh,DICpara.GaussPtOrder,beta,mu,USubpb1,FSubpb1,udual,vdual,alpha,mean(DICpara.winstepsize),0);
		[FSubpb2,~,~] = funGlobalNodalStrainQuadtree(DICmesh,USubpb2,DICpara.GaussPtOrder,0);
        ALSub2Time(ALSolveStep) = toc; toc
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ------- Smooth strain field --------
        if DICpara.DispSmoothness>1e-6, USubpb2 = funSmoothDispQuadtree(USubpb2,DICmesh,DICpara); end
        % ------- Don't change strain fields near the boundary --------
        for tempk=0:3, FSubpb2(4*DICmesh.markCoordHoleEdge-tempk) = FSubpb1(4*DICmesh.markCoordHoleEdge-tempk); end
        if DICpara.StrainSmoothness>1e-6, FSubpb2 = funSmoothStrainQuadtree(0.1*FSubpb2+0.9*FSubpb1,DICmesh,DICpara); end
        for tempk=0:1, USubpb2(2*markCoordHoleStrain-tempk) = USubpb1(2*markCoordHoleStrain-tempk); end
        for tempk=0:3, FSubpb2(4*markCoordHoleStrain-tempk) = FSubpb1(4*markCoordHoleStrain-tempk); end
         
		save(['Subpb2_step',num2str(ALSolveStep)],'USubpb2','FSubpb2');
         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute norm of UpdateY
        USubpb2_Old = load(['Subpb2_step',num2str(ALSolveStep-1)],'USubpb2');
        USubpb2_New = load(['Subpb2_step',num2str(ALSolveStep)],'USubpb2');
        USubpb1_Old = load(['Subpb1_step',num2str(ALSolveStep-1)],'USubpb1');
        USubpb1_New = load(['Subpb1_step',num2str(ALSolveStep)],'USubpb1');
        if (mod(ImgSeqNum-2,DICpara.ImgSeqIncUnit) ~= 0 && (ImgSeqNum>2)) || (ImgSeqNum < DICpara.ImgSeqIncUnit)
            UpdateY = norm((USubpb2_Old.USubpb2 - USubpb2_New.USubpb2), 2)/sqrt(size(USubpb2_Old.USubpb2,1));
            try
                UpdateY2 = norm((USubpb1_Old.USubpb1 - USubpb1_New.USubpb1), 2)/sqrt(size(USubpb1_Old.USubpb1,1));
            catch
            end
        end
        try
            disp(['Update local step  = ',num2str(UpdateY2)]);
            disp(['Update global step = ',num2str(UpdateY)]);
        catch
        end
        fprintf('*********************************** \n \n');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Update dual variables------------------------------
        udual = FSubpb2 - FSubpb1; vdual = USubpb2 - USubpb1; 
		 
        save(['uvdual_step',num2str(ALSolveStep)],'udual','vdual');
        try
        if UpdateY < tol2 || UpdateY2 < tol2
            break
        end
        catch
        end
         
    end
    fprintf('------------ Section 6 Done ------------ \n \n')
 
    % Save data
    ResultDisp{ImgSeqNum-1}.U = full(USubpb2);
    ResultDisp{ImgSeqNum-1}.ALSub1BadPtNum = ALSub1BadPtNum;
    ResultDefGrad{ImgSeqNum-1}.F = full(FSubpb2);  

    
    % % Save data
    % ResultDisp{ImgSeqNum-1}.U = full(USubpb1);
    % ResultDisp{ImgSeqNum-1}.ALSub1BadPtNum = ALSub1BadPtNum;
    % ResultDefGrad{ImgSeqNum-1}.F = full(FSubpb1);
    

end    


%% ------ Plot ------
USubpb2World = USubpb2; USubpb2World(2:2:end) = -USubpb2(2:2:end); FSubpb2World = FSubpb2; 
close all; Plotdisp_show(USubpb2World,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4),DICpara,'EdgeColor');
Plotstrain_show(FSubpb2World,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4),DICpara,'EdgeColor');

%%%%%%%% These lines can be used to plot updated DIC subset window sizes %%%%%%%
% winsize_xy = [winsize_x_List, winsize_y_List]'; winsize_xy = winsize_xy(:);
% Plotdisp_show(winsize_xy,DICmesh.coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4),DICpara,'EdgeColor');
        
% ------ Save results ------
% Find img name and save all the results 
[~,imgname,imgext] = fileparts(file_name{1,end});
results_name = ['results_',imgname,'_ws',num2str(DICpara.winsize),'_st',num2str(DICpara.winstepsize),'.mat'];
save(results_name, 'file_name','DICpara','DICmesh','ResultDisp','ResultDefGrad','ResultFEMesh','ResultFEMeshEachFrame','ALSub1Time','ALSub2Time','ALSolveStep');


%% Section 7: Check convergence
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section is to check convergence of ADMM
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('------------ Section 7 Start ------------ \n')
% ====== Check convergence ======
fprintf('***** Check convergence ***** \n');
ALSolveStep1 = min(6,ALSolveStep);
disp('==== uhat^(k) - u^(k) ====');
for ALSolveStep = 1:ALSolveStep1
    USubpb2 = load(['Subpb2_step',num2str(ALSolveStep )],'USubpb2');
    USubpb1 = load(['Subpb1_step',num2str(ALSolveStep )],'USubpb1');
    UpdateY = norm((USubpb2.USubpb2 - USubpb1.USubpb1), 2)/sqrt(length(USubpb2.USubpb2));
    disp(num2str(UpdateY));
end
disp('==== Fhat^(k) - F^(k) ====');
for ALSolveStep = 1:ALSolveStep1
    FSubpb1 = load(['Subpb1_step',num2str(ALSolveStep )],'FSubpb1');
    FSubpb2 = load(['Subpb2_step',num2str(ALSolveStep )],'FSubpb2');
    UpdateF = norm((FSubpb1.FSubpb1 - FSubpb2.FSubpb2), 2)/sqrt(length(FSubpb1.FSubpb1));
    disp(num2str(UpdateF));
end
disp('==== uhat^(k) - uhat^(k-1) ====');
for ALSolveStep = 2:ALSolveStep1
    USubpb2_Old = load(['Subpb2_step',num2str(ALSolveStep-1)],'USubpb2');
    USubpb2_New = load(['Subpb2_step',num2str(ALSolveStep)],'USubpb2');
    UpdateY = norm((USubpb2_Old.USubpb2 - USubpb2_New.USubpb2), 2)/sqrt(length(USubpb2.USubpb2));
    disp(num2str(UpdateY));
end
disp('==== udual^(k) - udual^(k-1) ====');
for ALSolveStep = 2:ALSolveStep1
    uvdual_Old = load(['uvdual_step',num2str(ALSolveStep-1)],'udual');
    uvdual_New = load(['uvdual_step',num2str(ALSolveStep)],'udual');
    UpdateW = norm((uvdual_Old.udual - uvdual_New.udual), 2)/sqrt(length(uvdual_Old.udual));
    disp(num2str(UpdateW));
end
disp('==== vdual^(k) - vdual^(k-1) ====');
for ALSolveStep = 2:ALSolveStep1
    uvdual_Old = load(['uvdual_step',num2str(ALSolveStep-1)],'vdual');
    uvdual_New = load(['uvdual_step',num2str(ALSolveStep)],'vdual');
    Updatev = norm((uvdual_Old.vdual - uvdual_New.vdual), 2)/sqrt(length(uvdual_Old.vdual));
    disp(num2str(Updatev));
end
fprintf('------------ Section 7 Done ------------ \n \n')

% ------ Delete temp files ------
%%%%% Comment START %%%%%
% Uncomment these lines to delete temporary files
%   for tempi = 1:ALSolveStep
%       file_name_Subpb1 = ['Subpb1_step',num2str(tempi),'.mat'];
%       file_name_Subpb2 = ['Subpb2_step',num2str(tempi),'.mat'];
%       file_name_dual = ['uvdual_step',num2str(tempi),'.mat'];
%       delete(file_name_Subpb1); delete(file_name_Subpb2); delete(file_name_dual);
%   end
%%%%% Comment END %%%%%

% ------ clear temp variables ------
clear a ALSub1BadPtNum ALSub1Timetemp atemp b btemp cc ConvItPerEletemp hbar Hbar 
clear coordinatesFEMQuadtree elementsFEMQuadtree 


%% ====== Transform "incremental" displacement fields to "cumulative" displacement fields ======
tempx = ResultFEMeshEachFrame{1}.coordinatesFEM(:,1);
tempy = ResultFEMeshEachFrame{1}.coordinatesFEM(:,2);
coord = [tempx,tempy]; coordCurr = coord;  
hbar = waitbar(0,'Calculate cumulative disp from incremental disp');

for ImgSeqNum = 2 :  length(ImgNormalized)
     
    waitbar((ImgSeqNum-1)/(size(file_name,2)-1));
    tempx = ResultFEMeshEachFrame{ImgSeqNum-1}.coordinatesFEM(:,1);
    tempy = ResultFEMeshEachFrame{ImgSeqNum-1}.coordinatesFEM(:,2);
    
    tempu = ResultDisp{ImgSeqNum-1}.U(1:2:end);
    tempv = ResultDisp{ImgSeqNum-1}.U(2:2:end);
    
    op2_x = rbfcreate( [tempx,tempy]',[tempu]','RBFFunction', 'thinplate');
    rbfcheck_maxdiff = rbfcheck(op2_x); 
    if rbfcheck_maxdiff > 1e-3, disp('Please check rbf interpolation! Pause here.'); pause; end
    disp_x = rbfinterp([coordCurr(:,1),coordCurr(:,2)]', op2_x );
    
    op2_y = rbfcreate( [tempx,tempy]',[tempv]','RBFFunction', 'thinplate');
    rbfcheck_maxdiff = rbfcheck(op2_y); 
    if rbfcheck_maxdiff > 1e-3, disp('Please check rbf interpolation! Pause here.'); pause; end
    disp_y = rbfinterp([coordCurr(:,1),coordCurr(:,2)]', op2_y );
    
    coordCurr = coordCurr + [disp_x(:), disp_y(:)];
    U_cum = (coordCurr - coord)'; U_cum = U_cum(:); 
    ResultDisp{ImgSeqNum-1}.U_cum_store = U_cum; % Store cumulative displacement field
    
end

close(hbar);
%%%%%%% TODO %%%%%%%
Plotdisp_show( ResultDisp{2-1}.U_cum_store,ResultFEMeshEachFrame{1}.coordinatesFEM,ResultFEMeshEachFrame{1}.elementsFEM(:,1:4),DICpara,'EdgeColor');



%% Section 8: Compute strains
fprintf('------------ Section 8 Start ------------ \n')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section is to compute strain fields and plot disp and strain results
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------ Convert units from pixels to the physical world ------
DICpara.um2px = funParaInput('ConvertUnit');
% ------ Smooth displacements ------
DICpara.DoYouWantToSmoothOnceMore = 1; % No need to smooth disp fields
DICpara.smoothness = funParaInput('RegularizationSmoothness'); % Regularization to smooth strain fields           
% ------ Choose strain computation method ------
DICpara.MethodToComputeStrain = 2; % funParaInput('StrainMethodOp'); 
if DICpara.MethodToComputeStrain == 2 % Compute strain method II: Use Plane Fitting method
    prompt = 'What is your half window size (unit: px): ';
    Rad = input(prompt);      
end
% ------ Choose strain type (infinitesimal, Eulerian, Green-Lagrangian) ------
DICpara.StrainType = funParaInput('StrainType');
% ------ Choose image to plot (first only, second and next images) ------
if length(ImgNormalized)==2, DICpara.Image2PlotResults = funParaInput('Image2PlotResults');
else DICpara.Image2PlotResults = 1; % Plot over current, deformed image by default
end
% ------ Save fig format ------
DICpara.MethodToSaveFig = funParaInput('SaveFigFormat');
% ------ Choose overlay image transparency ------
DICpara.OrigDICImgTransparency = 1; 
if DICpara.MethodToSaveFig == 1  
    DICpara.OrigDICImgTransparency = funParaInput('OrigDICImgTransparency');         
end


%% ====== Start main part ======
% This section is to calculate strain fields based on the transformed
% cumulative displacements: [F] = [D][U]
% [F] = [..., F11_nodei, F21_nodei, F12_nodei, F22_nodei, ...]';
% [u] = [..., U1_nodei, U2_nodei, ...]';
% [D]: finite difference/finite element operator to compute first derivatives

% v = VideoWriter('video_cav_DIC.mp4');
% v.FrameRate = 10;
% open(v);

 
for ImgSeqNum =  [ 2 :  length(ImgNormalized) ] 
    
    close all; disp(['Current image frame #: ', num2str(ImgSeqNum),'/',num2str(length(ImgNormalized))]);
     
    %%%%% Load deformed image %%%%%%%
    gNormalizedMask = double( ImgMask{ImgSeqNum} ); % Load the mask file of current deformed frame
    gNormalized = ImgNormalized{ImgSeqNum} .* gNormalizedMask ; % Load current deformed frame 
    Dg = funImgGradient(gNormalized,gNormalized,gNormalizedMask); % Finite difference to compute image grayscale gradients;
    
    fNormalizedMask = double( ImgMask{1} ); % 
    DICpara.ImgRefMask = fNormalizedMask;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ====== Old version ======
    % fNormalizedNewIndex = ImgSeqNum-mod(ImgSeqNum-2,DICpara.ImgSeqIncUnit)-1;
    % if DICpara.ImgSeqIncUnit > 1
    %     FEMeshIndLast = floor(fNormalizedNewIndex/DICpara.ImgSeqIncUnit);
    % elseif DICpara.ImgSeqIncUnit == 1
    %     FEMeshIndLast = floor(fNormalizedNewIndex/DICpara.ImgSeqIncUnit)-1;
    % end
    % FEMeshInd = FEMeshIndLast + 1;
    % 
    % if FEMeshInd == 1
    %     USubpb2 = ResultDisp{ImgSeqNum-1}.U; %+ ResultDisp{10}.U + ResultDisp{20}.U;
    %     coordinatesFEM = ResultFEMesh{1}.coordinatesFEM; 
    %     elementsFEM = ResultFEMesh{1}.elementsFEM;
    %     if (ImgSeqNum-1 == 1) || (DICpara.ImgSeqIncROIUpdateOrNot==1), UFEMesh = 0*USubpb2; end
    % else
    %     USubpb2 = ResultDisp{ImgSeqNum-1}.U;
    %     if mod(ImgSeqNum-2,DICpara.ImgSeqIncUnit) == 0
    %         coordinatesFEM = ResultFEMesh{FEMeshInd}.coordinatesFEM;
    %         elementsFEM = ResultFEMesh{FEMeshInd}.elementsFEM;
    %         coordinatesFEMLast = ResultFEMesh{FEMeshIndLast}.coordinatesFEM;
    %         UFEMeshLast = ResultDisp{ImgSeqNum-2}.U + UFEMesh;
    %         xq = coordinatesFEM(:,1); yq = coordinatesFEM(:,2);
    %         UFEMesh = 0*USubpb2;
    %         UFEMesh(1:2:end) = griddata(coordinatesFEMLast(:,1),coordinatesFEMLast(:,2),UFEMeshLast(1:2:end),xq,yq,'v4');
    %         UFEMesh(2:2:end) = griddata(coordinatesFEMLast(:,1),coordinatesFEMLast(:,2),UFEMeshLast(2:2:end),xq,yq,'v4');
    %     end
    %     USubpb2 = USubpb2 + UFEMesh;
    % end
    %
    % FSubpb2 = ResultDefGrad{ImgSeqNum-1}.F;
    % coordinatesFEM = ResultFEMeshEachFrame{ImgSeqNum-1}.coordinatesFEM;
    % elementsFEM = ResultFEMeshEachFrame{ImgSeqNum-1}.elementsFEM;
    
    USubpb2 = ResultDisp{ImgSeqNum-1}.U_cum_store;
    coordinatesFEM = ResultFEMeshEachFrame{1}.coordinatesFEM;
    elementsFEM = ResultFEMeshEachFrame{1}.elementsFEM;
 
    try markCoordHoleEdge = ResultFEMeshEachFrame{ImgSeqNum-1}.markCoordHoleEdge; catch; end
    DICmesh.coordinatesFEM = coordinatesFEM;
    DICmesh.elementsFEM = elementsFEM;
    coordinatesFEMWorld = DICpara.um2px*[coordinatesFEM(:,1),size(ImgNormalized{1},2)+1-coordinatesFEM(:,2)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    % ------ Plotting and Compute Strain-------
    if size(USubpb2,1) == 1
        ULocal = USubpb2_New.USubpb2; % FLocal = FSubpb2.FSubpb2; 
    else
        ULocal = USubpb2; % FLocal = FSubpb2;
    end
    UWorld = DICpara.um2px*ULocal; UWorld(2:2:end) = -UWorld(2:2:end);   % close all; Plotuv(UWorld,x0,y0World);
     
    % ------ Smooth displacements ------
    %prompt = 'Do you want to smooth displacement? (0-yes; 1-no)';
    %DoYouWantToSmoothOnceMore = input(prompt); 
    SmoothTimes = 0;
    try
        while DICpara.DoYouWantToSmoothOnceMore == 0 && SmoothTimes < 3
            ULocal = funSmoothDispQuadtree(ULocal,DICmesh,DICpara);
            %close all; Plotuv(ULocal,x0,y0); %DICpara.DoYouWantToSmoothOnceMore = input(prompt);
            SmoothTimes = SmoothTimes + 1;
        end
    catch
    end
    
    % ----- Compute strain field ------
    ComputeStrainQuadtree;  
    
    % ------ Smooth strain fields ------
%     DICpara.DoYouWantToSmoothOnceMore = 0; % Smooth strain fields if necessary
%     %prompt = 'Do you want to smooth displacement? (0-yes; 1-no)';
%     %DoYouWantToSmoothOnceMore = input(prompt); 
%     SmoothTimes = 0;
%     try
%         while DICpara.DoYouWantToSmoothOnceMore == 0 && SmoothTimes < 3
%             FStraintemp = funSmoothStrainQuadtree(FStraintemp,DICmesh,DICpara);
%             %close all; Plotuv(ULocal,x0,y0); %DICpara.DoYouWantToSmoothOnceMore = input(prompt);
%             SmoothTimes = SmoothTimes + 1;
%         end
%     catch
%     end
%     FStrainWorld = FStraintemp; FStrainWorld(2:4:end) = -FStrainWorld(2:4:end); FStrainWorld(3:4:end) = -FStrainWorld(3:4:end); 
     
     
    % ------ Plot disp and strain ------
    if DICpara.OrigDICImgTransparency == 1
        Plotdisp_show(UWorld,coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4),DICpara,'NoEdgeColor');
        [strain_exx,strain_exy,strain_eyy,strain_principal_max,strain_principal_min,strain_maxshear,strain_vonMises] = ...
                   Plotstrain0Quadtree(FStrainWorld,coordinatesFEMWorld,elementsFEM(:,1:4),DICpara);
    
    else % Plot over raw DIC images
        if DICpara.Image2PlotResults == 0 % Plot over the first image; "file_name{1,1}" corresponds to the first image	
            PlotdispQuadtree(UWorld,coordinatesFEMWorld,elementsFEM(:,1:4),file_name{1,1},DICpara);
            [strain_exx,strain_exy,strain_eyy,strain_principal_max,strain_principal_min, ...
                strain_maxshear,strain_vonMises] = PlotstrainQuadtree(UWorld,FStrainWorld, ...
                coordinatesFEMWorld,elementsFEM(:,1:4),file_name{1,1},DICpara);

        else % Plot over second or next deformed images
            
            %%%%%% Old codes: without applying mask files %%%%%%
            % PlotdispQuadtree(UWorld,coordinatesFEMWorld,elementsFEM(:,1:4),...
            %    file_name{1,ImgSeqNum},DICpara);
            % 
            % [strain_exx,strain_exy,strain_eyy,strain_principal_max,strain_principal_min, ...
            %    strain_maxshear,strain_vonMises] = PlotstrainQuadtree(UWorld,FStraintemp, ...
            %    coordinatesFEMWorld,elementsFEM(:,1:4),file_name{1,ImgSeqNum},DICpara);

            %%%%%% New codes: applying mask files %%%%%%
            PlotdispQuadtreeMasks(UWorld,coordinatesFEMWorld,elementsFEM(:,1:4),...
                file_name{1, ImgSeqNum}, ...
                ImgMask{ ImgSeqNum },DICpara);
            
            [strain_exx,strain_exy,strain_eyy,strain_principal_max,strain_principal_min, ...
               strain_maxshear,strain_vonMises] = PlotstrainQuadtreeMasks(UWorld,FStraintemp, ...
               coordinatesFEMWorld,elementsFEM(:,1:4),file_name{1, ImgSeqNum}, ...
               ImgMask{ ImgSeqNum },DICpara);
           

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO 
%             PlotdispQuadtreePolarMasks(UWorld,coordinatesFEMWorld,elementsFEM(:,1:4),...
%                 file_name{1, ImgSeqNum}, ...
%                 ImgMask{ ImgSeqNum },DICpara,ImgSeqNum);
            
%             [strain_err,strain_ert,strain_ett,strain_principal_max,strain_principal_min, ...
%                strain_maxshear,strain_vonMises] = PlotstrainQuadtreePolarMasks(UWorld,FStrainWorld, ...
%                coordinatesFEMWorld,elementsFEM(:,1:4),file_name{1, ImgSeqNum}, ...
%                ImgMask{ ImgSeqNum },DICpara, [] ,ImgSeqNum);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TODO           
           

 %%%%%%%%% cav-dic part %%%%%%%%%% TODO 
 %            disp_u = UWorld(1:2:end); disp_v = UWorld(2:2:end);
%            coordinatesFEMWorldDef = [coordinatesFEMWorld(:,1)+ disp_u, ...
%                                      coordinatesFEMWorld(:,2)+ disp_v];
%  

%          bubble_y = (251-mean(CircleFitPar(end-20:end-1,1)))*DICpara.um2px;
%          bubble_x = mean(CircleFitPar(end-20:end-1,2))*DICpara.um2px;
%          
%          r = sqrt( (coordinatesFEMWorldDef(:,1)-bubble_x).^2 + (coordinatesFEMWorldDef(:,2)-bubble_y).^2  );
%          theta = atan2(  (coordinatesFEMWorldDef(:,2)-bubble_y), coordinatesFEMWorldDef(:,1)-bubble_x);
%          
%          disp_r = cos(theta).*disp_u + sin(theta).*disp_v; % JY!!! correction on 08/15/2021
%          disp_t = - sin(theta).*disp_u + cos(theta).*disp_v;
%          
%       
%          strain_logEtt = log(1+strain_ett);
%          strain_logErr = log(1+strain_err);
%          
%          Jacobian = (1+strain_err).*(1+strain_ett).^2; 
%          
%          frame = getframe(gcf);
%          writeVideo(v,frame);
%     
%         %%%%%%%%% cav-dic part %%%%%%%%%%

        end
    end
  
    % ----- Save strain results ------
    ResultStrain{ImgSeqNum-1} = struct('strainxCoord',coordinatesFEMWorld(:,1),'strainyCoord',coordinatesFEMWorld(:,2), ...
            'dispu',UWorld(1:2:end),'dispv',UWorld(2:2:end), ...
            'dudx',FStraintemp(1:4:end),'dvdx',FStraintemp(2:4:end),'dudy',FStraintemp(3:4:end),'dvdy',FStraintemp(4:4:end), ...
            'strain_exx',strain_exx,'strain_exy',strain_exy,'strain_eyy',strain_eyy, ...
            'strain_principal_max',strain_principal_max,'strain_principal_min',strain_principal_min, ...
            'strain_maxshear',strain_maxshear,'strain_vonMises',strain_vonMises);
     
%     ResultCavDIC{ImgSeqNum-1} = struct('bubble_center_x',bubble_x,'bubble_center_y',bubble_y, ...
%             'r',r,'theta',theta,'disp_r',disp_r,'disp_t',disp_t, ...
%             'R',Rnew(ImgSeqNum)*DICpara.um2px,'DICwinsizePhy',DICpara.winsize*DICpara.um2px,  ...
%             'dudx',dudx,'dvdx',dvdx,'dudy',dudy,'dvdy',dvdy, ...
%             'strain_err',strain_err,'strain_ert',strain_ert,'strain_ett',strain_ett, ...
%             'strain_logErr',strain_logErr,'strain_logEtt',strain_logEtt, ...
%             'Jacobian',Jacobian');


% % ------ Save figures for tracked displacement and strain fields ------
[~,imgname,imgext] = fileparts(file_name{1,ImgSeqNum}); % Find img name
  SaveFigFilesDispAndStrainQuadtree;
   % pause;
    
    
end
% ------ END of for-loop {ImgSeqNum = 2:length(ImgNormalized)} ------
fprintf('------------ Section 8 Done ------------ \n \n')


% ------ Save data again including solved strain fields ------
results_name = ['results_',imgname,'_ws',num2str(DICpara.winsize),'_st',num2str(DICpara.winstepsize),'.mat'];
save(results_name, 'file_name','DICpara','DICmesh','ResultDisp','ResultDefGrad','ResultFEMesh','ResultFEMeshEachFrame',...
                   'ALSub1Time','ALSub2Time','ALSolveStep','ResultStrain');

% save('results_cav_dic.mat', 'ResultCavDIC');
% close(v);


%% Section 9: Compute stress
fprintf('------------ Section 9 Start ------------ \n')
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This section is to compute stress fields and plot stress fields
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------ Choose material model ------ 
DICpara.MaterialModel = funParaInput('MaterialModel');
% ------ Define parameters in material models ------
if (DICpara.MaterialModel == 1) || (DICpara.MaterialModel == 2) % Linear elasticity
    fprintf('Define Linear elasticity parameters \n')
    fprintf("Young's modulus (unit: Pa). \n "); prompt = 'Input here (e.g., 69e9): '; 
    DICpara.MaterialModelPara.YoungsModulus = input(prompt); 
    fprintf("Poisson's ratio \n"); prompt = 'Input here (e.g., 0.3): '; 
    DICpara.MaterialModelPara.PoissonsRatio = input(prompt);
    fprintf('------------------------------------- \n');
end

% ------ Start main part ------
for ImgSeqNum = 2 : length(ImgNormalized)
    
    disp(['Current image frame #: ', num2str(ImgSeqNum),'/',num2str(length(ImgNormalized))]); close all;
    
    coordinatesFEM = ResultFEMeshEachFrame{ImgSeqNum-1}.coordinatesFEM;
    elementsFEM = ResultFEMeshEachFrame{ImgSeqNum-1}.elementsFEM;
    coordinatesFEMWorldDef = DICpara.um2px*[coordinatesFEM(:,1),size(ImgNormalized{1},2)+1-coordinatesFEM(:,2)] + ...
                             DICpara.Image2PlotResults*[ResultStrain{ImgSeqNum-1}.dispu, ResultStrain{ImgSeqNum-1}.dispv];
     
    % ------ Plot stress ------
    if DICpara.OrigDICImgTransparency == 1
        [stress_sxx,stress_sxy,stress_syy, stress_principal_max_xyplane, ...
                stress_principal_min_xyplane, stress_maxshear_xyplane, ...
                stress_maxshear_xyz3d, stress_vonMises]  =  Plotstress0Quadtree( ...
            DICpara,ResultStrain{ImgSeqNum-1},coordinatesFEMWorldDef,elementsFEM(:,1:4)); 
        
    else % Plot over raw DIC images
        if DICpara.Image2PlotResults == 0 % Plot over the first image; "file_name{1,1}" corresponds to the first image	
            [stress_sxx,stress_sxy,stress_syy, stress_principal_max_xyplane, ...
                stress_principal_min_xyplane, stress_maxshear_xyplane, ...
                stress_maxshear_xyz3d, stress_vonMises] = PlotstressQuadtree( ...
                DICpara,ResultStrain{ImgSeqNum-1},coordinatesFEMWorldDef,elementsFEM(:,1:4),file_name{1,1});
             
        else % Plot over second or next deformed images
           [stress_sxx,stress_sxy,stress_syy, stress_principal_max_xyplane, ...
                stress_principal_min_xyplane, stress_maxshear_xyplane, ...
                stress_maxshear_xyz3d, stress_vonMises] = PlotstressQuadtree( ...
                DICpara,ResultStrain{ImgSeqNum-1},coordinatesFEMWorldDef,elementsFEM(:,1:4),file_name{1,ImgSeqNum});
 
        end
    end
    
    
    % ------ Save figures for computed stress fields ------
    SaveFigFilesStress; 
    
    % ----- Save strain results ------
    ResultStress{ImgSeqNum-1} = struct('stressxCoord',ResultStrain{ImgSeqNum-1}.strainxCoord,'stressyCoord',ResultStrain{ImgSeqNum-1}.strainyCoord, ...
            'stress_sxx',stress_sxx,'stress_sxy',stress_sxy,'stress_syy',stress_syy, ...
            'stress_principal_max_xyplane',stress_principal_max_xyplane, 'stress_principal_min_xyplane',stress_principal_min_xyplane, ...
            'stress_maxshear_xyplane',stress_maxshear_xyplane,'stress_maxshear_xyz3d',stress_maxshear_xyz3d, ...
            'stress_vonMises',stress_vonMises);
        
end
% ------ END of for-loop {ImgSeqNum = 2:length(ImgNormalized)} ------
fprintf('------------ Section 9 Done ------------ \n \n')

% ------ Save data again including solved stress fields ------
results_name = ['results_',imgname,'_ws',num2str(DICpara.winsize),'_st',num2str(DICpara.winstepsize),'.mat'];
save(results_name, 'file_name','DICpara','DICmesh','ResultDisp','ResultDefGrad','ResultFEMesh','ResultFEMeshEachFrame', ...
                   'ALSub1Time','ALSub2Time','ALSolveStep','ResultStrain','ResultStress');


%% Section 10: Plot the generated quadtree mesh 
v = VideoWriter('video_mesh.mp4','MPEG-4');
v.FrameRate = 5;
open(v);
figure,
for ImgSeqNum = 2 : (1+size(ResultDisp,1))
    
    clf; patch('Faces', DICmesh.elementsFEM(:,1:4), 'Vertices', DICmesh.coordinatesFEMWorld + ...
        [ResultDisp{ImgSeqNum-1}.U(1:2:end), -ResultDisp{ImgSeqNum-1}.U(2:2:end)], 'Facecolor','none','linewidth',1)
    xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
    tt = title(['Frame #',num2str(ImgSeqNum)],'fontweight','normal');
    set(tt,'Interpreter','latex','fontsize',10);
    axis equal; axis tight; set(gca,'fontsize',18); set(gcf,'color','w'); box on;
    a = gca; a.TickLabelInterpreter = 'latex';
    
    frame = getframe(gcf);
    writeVideo(v,frame);
    
end
close(v);


%% For Nomin's dataset
exx_cum_store = []; exx_std_cum_store = [];
exy_cum_store = []; exy_std_cum_store = [];
eyy_cum_store = []; eyy_std_cum_store = [];

for ImgSeqNum = 2:length(file_name)
    
    dudx = ResultStrain{ImgSeqNum-1}.dudx;
    dvdx = ResultStrain{ImgSeqNum-1}.dvdx;
    dudy = ResultStrain{ImgSeqNum-1}.dudy;
    dvdy = ResultStrain{ImgSeqNum-1}.dvdy;
    
    exx_cum_store(ImgSeqNum-1) = mean(dudx);
    exy_cum_store(ImgSeqNum-1) = mean(0.5*(dudy+dvdx));
    eyy_cum_store(ImgSeqNum-1) = mean(dvdy);
    
    exx_std_cum_store(ImgSeqNum-1) = std(dudx);
    exy_std_cum_store(ImgSeqNum-1) = std(0.5*(dudy+dvdx));
    eyy_std_cum_store(ImgSeqNum-1) = std(dvdy);
    
end

figure, errorbar([2:length(file_name)],exx_cum_store,exx_std_cum_store,'s-');
hold on;  errorbar([2:length(file_name)],exy_cum_store,exy_std_cum_store,'^-');
hold on;  errorbar([2:length(file_name)],eyy_cum_store,eyy_std_cum_store,'o-');


set(gca,'fontsize',20);
xlabel('Frame #');
ylabel('Strain');
legend('exx','exy','eyy');


