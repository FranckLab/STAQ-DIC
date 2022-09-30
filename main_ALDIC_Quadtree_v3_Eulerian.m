% ---------------------------------------------
% Augmented Lagrangian Digital Image Correlation (ALDIC_Quadtree_v2)
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
addpath("./func",'./plotFiles','./func_quadtree','./func_quadtree/refinement','./plotFiles/export_fig-d966721/'); 
addpath('./Images_Quadtree_demo/Images_Sample12'); % TODO: addpath("./YOUR IMAGE FOLDER"); 
addpath('./rbfinterp');
fprintf('------------ Section 1 Done ------------ \n \n')


%% Section 2: Load DIC parameters and set up DIC parameters 
fprintf('------------ Section 2 Start ------------ \n')
% ====== Read images ====== 
close all; [file_name,Img,DICpara] = ReadImageQuadtree; % Load DIC raw images

% ====== Read image mask files ======
[mask_file_name,ImgMask] = ReadImageMasks; % Load DIC image mask files

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
fprintf('------------ Section 2 Done ------------ \n \n')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To solve each frame in an image sequence
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ImgSeqNum = 2 : length(ImgNormalized) 
 

    disp(['Current image frame #: ', num2str(ImgSeqNum),'/',num2str(length(ImgNormalized))]);
    
    % ====== Load image mask files & Compute image gradients ======
    %%%%% Load reference image %%%%%
    fNormalizedMask = double( ImgMask{1} ) ;
    fNormalized = ImgNormalized{1}; % .* fNormalizedMask; % Load the first referece image 
    % Df = funImgGradient(fNormalized,fNormalized,fNormalizedMask); % Finite difference to compute image grayscale gradients;
    
    %%%%% Load deformed image %%%%%%%
    gNormalized = ImgNormalized{ ImgSeqNum } ; % Load current deformed image frame 
    gNormalizedMask = double( ImgMask{ImgSeqNum} );
    % gNormalized = gNormalized .* gNormalizedMask ;
    Dg = funImgGradient(gNormalized,gNormalized,gNormalizedMask); % Finite difference to compute image grayscale gradients;
    
    DICpara.ImgRefMask = gNormalizedMask;
    
    figure, subplot(2,2,1); imshow(fNormalized'); title('Ref img f')
    subplot(2,2,2); imshow(gNormalized'); title('Def img g')
    subplot(2,2,3); imshow(fNormalizedMask'); title('f mask')
    subplot(2,2,4); imshow(gNormalizedMask'); title('g mask')
    
    
    %% Section 3: Compute an initial guess of the unknown displacement field
    fprintf('\n'); fprintf('------------ Section 3 Start ------------ \n')
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This section is to find or update an initial guess of the unknown displacements.
    % The key idea is to either to use a new FFT-based cross correlation peak fitting,
    % or use the results from previous frames to compute a new initial guess for the next frame;
    % Particularly in the incremental mode DIC, the reference image can also be updated, e.g.,
    % " fNormalized = ImgNormalized{ImgSeqNum-mod(ImgSeqNum-1,ImgSeqIncUnit)}; "
    %
    % DICpara.NewFFTSearch = 0; % If you want to apply the FFT-based cross correlation to 
    % compute the initial guess for each frame, please make sure that "DICpara.NewFFTSearch = 0". 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%% One practical strategy is to let first 7 frames do the FFT-based 
    %%%%% cross correlation and then let data driven method to estimate new 
    %%%%% initial guesses for other frames 
    if ImgSeqNum < 7
        DICpara.NewFFTSearch = 1; % Use FFT-based cross correlation to compute the initial guess
    else
        DICpara.NewFFTSearch = 0; % Apply data driven method to estimate initial guesses for later frames
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ImgSeqNum == 2 || DICpara.NewFFTSearch == 1 % Apply FFT-based cross correlation to compute the initial guess 
         
        % ====== FFT-based cross correlation ======
        [DICpara,x0temp,y0temp,u,v,cc]= IntegerSearch(gNormalized,fNormalized,file_name,DICpara);
        % ====== DIC uniform FE-mesh set up ======
        [DICmesh] = MeshSetUp(x0temp,y0temp,DICpara); clear x0temp y0temp;
        % ====== Initial Value ======
        U0 = Init(u,v,cc.max,DICmesh.x0,DICmesh.y0,0); % [Temp code:] PlotuvInit;
         
        % ====== Deal with incremental mode (i) ======
        fNormalizedNewIndex = ImgSeqNum-mod(ImgSeqNum-2,DICpara.ImgSeqIncUnit)-1;
        if DICpara.ImgSeqIncUnit == 1, fNormalizedNewIndex = fNormalizedNewIndex-1; end
         
        % ====== Generate a quadtree mesh considering sample's complex geometry ======
        DICmesh.elementMinSize = 4; % min element size in the refined quadtree mesh
        GenerateQuadtreeMeshEulerian; % Generate a quadtree mesh
        DICmesh.coordinatesFEM_Eulerian = DICmesh.coordinatesFEM; 
        DICmesh.coordinatesFEM_EulerianWorld = [DICmesh.coordinatesFEM_Eulerian(:,1),size(DICpara.ImgRefMask,2)+1-DICmesh.coordinatesFEM_Eulerian(:,2)];
        
        % ====== Deal with incremental mode (ii) ======
        ResultFEMesh{1+floor(fNormalizedNewIndex/DICpara.ImgSeqIncUnit)} = ... % To save first mesh info
            struct( 'coordinatesFEM_Eulerian',DICmesh.coordinatesFEM_Eulerian,'elementsFEM',DICmesh.elementsFEM, ...
            'winsize',DICpara.winsize,'winstepsize',DICpara.winstepsize,'gridxyROIRange',DICpara.gridxyROIRange );
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif mod(ImgSeqNum-2,DICpara.ImgSeqIncUnit) == 0 % To update ref image in incremental mode
        fNormalizedNewIndex = ImgSeqNum-mod(ImgSeqNum-2,DICpara.ImgSeqIncUnit)-1;
        if DICpara.ImgSeqIncUnit == 1,  fNormalizedNewIndex = fNormalizedNewIndex-1; end
        fNormalized = ImgNormalized{fNormalizedNewIndex}; % Update reference
        [DICpara,DICmesh] = ReadImageRefUpdate(file_name,ImgSeqNum,ResultDisp{ImgSeqNum-2}.U,DICpara,DICmesh); % Update reference image if needed;
        U0 = zeros(2*size(DICmesh.coordinatesFEM_Eulerian,1),1); % [Temporary code: " PlotuvInit; "] 
        ResultFEMesh{1+floor(fNormalizedNewIndex/DICpara.ImgSeqIncUnit)} = ... % To save first mesh info
            struct( 'coordinatesFEM_Eulerian',DICmesh.coordinatesFEM_Eulerian,'elementsFEM',DICmesh.elementsFEM, ...
            'winsize',DICpara.winsize,'winstepsize',DICpara.winstepsize,'gridxyROIRange',DICpara.gridxyROIRange );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else % Use the solved results from the last frame as the new initial guess
        if ImgSeqNum < 7 % Import previous U for ImgSeqNum [2,6] 
            U0 = ResultDisp{ImgSeqNum-2}.U;
             
        else % When ImgSeqNum > 6: POD predicts next disp U0 from previous results of (ImgSeqNum+[-5:1:-1])
            nTime = 5; 
            
            xnodes = DICpara.winstepsize+DICpara.gridxyROIRange.gridx(1) : DICpara.winstepsize : -DICpara.winstepsize+DICpara.gridxyROIRange.gridx(2);
            ynodes = DICpara.winstepsize+DICpara.gridxyROIRange.gridy(1) : DICpara.winstepsize : -DICpara.winstepsize+DICpara.gridxyROIRange.gridy(2);
            [x0temp,y0temp] = ndgrid(xnodes,ynodes);
            np = size(x0temp,1)*size(x0temp,2);
            T_data_u = zeros(nTime,np); T_data_v = zeros(nTime,np);
            
            for tempi =  1:nTime
                
                tempx = ResultFEMeshEachFrame{ImgSeqNum-(2+nTime)+tempi, 1}.coordinatesFEM_Eulerian(:,1);
                tempy = ResultFEMeshEachFrame{ImgSeqNum-(2+nTime)+tempi, 1}.coordinatesFEM_Eulerian(:,2);
                tempu = -ResultDisp{ImgSeqNum-(2+nTime)+tempi, 1}.U(1:2:end);
                tempv = -ResultDisp{ImgSeqNum-(2+nTime)+tempi, 1}.U(2:2:end);
                 
                op1 = rbfcreate( [tempx(:),tempy(:)]',[tempu(:)]','RBFFunction','thinplate','RBFSmooth',1); rbfcheck(op1);
                unodes = rbfinterp([x0temp(:),y0temp(:)]', op1 );
                op2 = rbfcreate( [tempx(:),tempy(:)]',[tempv(:)]','RBFFunction','thinplate','RBFSmooth',1); rbfcheck(op2);
                vnodes = rbfinterp([x0temp(:),y0temp(:)]', op2 );       
  
                T_data_u(tempi,:) = unodes(:);
                T_data_v(tempi,:) = vnodes(:);
                
            end
            nB = 3; t_train = [ImgSeqNum-1-nTime:ImgSeqNum-2]'; t_pre = [ImgSeqNum-1]';
            [u_pred,~,~,~] = funPOR_GPR(T_data_u,t_train,t_pre,nB);
            [v_pred,~,~,~] = funPOR_GPR(T_data_v,t_train,t_pre,nB);
            trained_u = reshape(u_pred(1,:),size(x0temp)); 
            trained_v = reshape(v_pred(1,:),size(x0temp));
            U0 = [trained_u(:), trained_v(:)]'; U0 = U0(:);
            
            
            % ====== DIC uniform FE-mesh set up ======
            [DICmesh] = MeshSetUp(x0temp,y0temp,DICpara); % clear x0temp y0temp;
            % ====== Initial Value ======
            U0 = Init(trained_u,trained_v,[],DICmesh.x0,DICmesh.y0,0); % [Temp code:] PlotuvInit;
        
            %%%%%% U0 = [trained_u(:),trained_v(:)]'; U0 = U0(:);
            
            for tempi = 1:size(trained_u,1)
                for tempj = 1:size(trained_u,2)
                    try
                        if ~gNormalizedMask(x0temp(tempi,tempj),y0temp(tempi,tempj))  
                            
                            U0(2*(tempj+(tempi-1)*(size(u,2)))) = nan;
                            U0(2*(tempj+(tempi-1)*(size(u,2)))-1) = nan;
                            
                        end
                    catch
                    end
                end
            end
            
            % ====== Generate a quadtree mesh considering sample's complex geometry ======
            DICmesh.elementMinSize = 4; % min element size in the refined quadtree mesh
            GenerateQuadtreeMeshEulerian; % Generate a quadtree mesh
            DICmesh.coordinatesFEM_Eulerian = DICmesh.coordinatesFEM; 
            DICmesh.coordinatesFEM_EulerianWorld = [DICmesh.coordinatesFEM_Eulerian(:,1),size(DICpara.ImgRefMask,2)+1-DICmesh.coordinatesFEM_Eulerian(:,2)];
        
        end
    end
    
    % ====== Compute f(X)-g(x+u) ======
    % PlotImgDiff(x0,y0,u,v,fNormalized,gNormalized);
    fprintf('------------ Section 3 Done ------------ \n \n')

 
    %% Section 4: ALDIC Subproblem 1 -or- Local ICGN Subset DIC
    fprintf('------------ Section 4 Start ------------ \n')
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This section is to solve the first local step in ALDIC: Subproblem 1
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    % ====== ALStep 1 Subproblem1: Local Subset DIC ======
    mu=0; beta=0; tol=1e-2; ALSolveStep=1; ALSub1Time=zeros(6,1); ALSub2Time=zeros(6,1); 
    ConvItPerEle=zeros(size(DICmesh.coordinatesFEM_Eulerian,1),6); ALSub1BadPtNum=zeros(6,1);
    disp(['***** Start step',num2str(ALSolveStep),' Subproblem1 *****'])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ------ Start Local DIC IC-GN iteration ------
    [USubpb1_Eulerian,FSubpb1_Eulerian,HtempPar,ALSub1Timetemp,ConvItPerEletemp,LocalICGNBadPtNumtemp] = ...
        LocalICGN(U0,DICmesh.coordinatesFEM_Eulerian,Dg,gNormalized,fNormalized,DICpara,'GaussNewton',tol);
    ALSub1Time(ALSolveStep) = ALSub1Timetemp; ConvItPerEle(:,ALSolveStep) = ConvItPerEletemp; ALSub1BadPtNum(ALSolveStep) = LocalICGNBadPtNumtemp; toc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ------  Manually find some bad points from Local Subset ICGN step ------
    % Comment these lines below if you don't have local bad points
    % %%%%% Comment START %%%%%%
    %  [USubpb1_Eulerian,FSubpb1_Eulerian] = funRemoveOutliersQuadtree(DICmesh,DICpara,USubpb1_Eulerian,FSubpb1_Eulerian);
    %  disp('--- Remove bad points done ---')
    % %%%%% Comment END %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    % ------ Plot ------
    USubpb1_EulerianWorld = USubpb1_Eulerian; USubpb1_EulerianWorld(2:2:end) = -USubpb1_Eulerian(2:2:end); FSubpb1_EulerianWorld = FSubpb1_Eulerian; 
    close all; Plotdisp_show(USubpb1_EulerianWorld,DICmesh.coordinatesFEM_EulerianWorld,DICmesh.elementsFEM(:,1:4),DICpara,'EdgeColor');
    Plotstrain_show(FSubpb1_EulerianWorld,DICmesh.coordinatesFEM_EulerianWorld,DICmesh.elementsFEM(:,1:4),DICpara,'EdgeColor');
    save(['Subpb1_step',num2str(ALSolveStep)],'USubpb1_Eulerian','FSubpb1_Eulerian');
    fprintf('------------ Section 4 Done ------------ \n \n')

    %%%%% Plot ICGN iteration step %%%%%
    % Convtemp = [ConvItPerEletemp]'; Convtemp = Convtemp(:);
    % PlotMesh_show(Convtemp,DICmesh.coordinatesFEM_EulerianWorld,DICmesh.elementsFEM(:,1:4),DICpara,'EdgeColor');
    
   
    %% Section 5: Subproblem 2 -- solve the global compatible displacement field
    fprintf('------------ Section 5 Start ------------ \n'); tic;
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This section is to solve the global step in ALDIC Subproblem 2
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    % ======= ALStep 1 Subproblem 2: Global constraint =======
    % ------ Smooth displacements for a better F ------
    DICpara.DispFilterSize=0; DICpara.DispFilterStd=0; DICpara.StrainFilterSize=0; DICpara.StrainFilterStd=0; LevelNo=1;
    DICpara.DispSmoothness = 0; DICpara.StrainSmoothness = 1;
    if DICpara.DispSmoothness>1e-6, USubpb1_Eulerian = funSmoothDispQuadtreeRBF(USubpb1_Eulerian,DICmesh,DICpara); end
    if DICpara.StrainSmoothness>1e-6, FSubpb1_Eulerian = funSmoothStrainQuadtreeRBF(FSubpb1_Eulerian,DICmesh,DICpara); end
    
    % close all; Plotdisp_show(USubpb1_Eulerian2,DICmesh.coordinatesFEM_EulerianWorld,DICmesh.elementsFEM(:,1:4),DICpara,'EdgeColor');
    % Plotstrain_show(FSubpb2_Eulerian ,DICmesh.coordinatesFEM_EulerianWorld,DICmesh.elementsFEM(:,1:4),DICpara,'EdgeColor');
    
    
	% ====== Define penalty parameter ======
    mu = 1e-3; udual = 0*FSubpb1_Eulerian; vdual = 0*USubpb1_Eulerian; 
    betaList = [1e-3,1e-2,1e-1]*mean(DICpara.winstepsize).^2.*mu; % Tune beta in the betaList 
    Err1 = zeros(length(betaList),1); Err2 = Err1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp(['***** Start step',num2str(ALSolveStep),' Subproblem2 *****']);
    DICpara.GaussPtOrder = 2; alpha = 0;  % No regularization added
    % ====== Solver using finite element method ======
    if ImgSeqNum == 2
        for tempk = 1:length(betaList)
           
            beta = betaList(tempk); display(['Try #',num2str(tempk),' beta = ',num2str(beta)]);
            GaussPtOrder=3; alpha=0; [USubpb2_Eulerian] = Subpb2Quadtree(DICmesh,DICpara.GaussPtOrder,beta,mu,USubpb1_Eulerian,FSubpb1_Eulerian,udual,vdual,alpha,mean(DICpara.winstepsize),0);
            % [FSubpb2_Eulerian,~,~] = funGlobalNodalStrainQuadtree(DICmesh,USubpb2_Eulerian,DICpara.GaussPtOrder,0);
            [FSubpb2_Eulerian] = funGlobalNodalStrainRBF(DICmesh,DICpara,USubpb2_Eulerian);
            % [FSubpb2_Eulerian] = funGlobalNodalStrainRBF(DICmesh,DICpara,USubpb1_Eulerian);
            % Plotstrain_show(FSubpb2_Eulerian ,DICmesh.coordinatesFEM_Eulerian ,DICmesh.elementsFEM(:,1:4),DICpara,'EdgeColor');
     
            Err1(tempk) = norm(USubpb1_Eulerian-USubpb2_Eulerian,2);
            Err2(tempk) = norm(FSubpb1_Eulerian-FSubpb2_Eulerian,2);
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
        [USubpb2_Eulerian] = Subpb2Quadtree(DICmesh,DICpara.GaussPtOrder,beta,mu,USubpb1_Eulerian,FSubpb1_Eulerian,udual,vdual,alpha,mean(DICpara.winstepsize),0);
        %[FSubpb2_Eulerian,~,~] = funGlobalNodalStrainQuadtree(DICmesh,USubpb2_Eulerian,DICpara.GaussPtOrder,0);
        [FSubpb2_Eulerian] = funGlobalNodalStrainRBF(DICmesh,DICpara,USubpb2_Eulerian);
        ALSub2Time(ALSolveStep) = toc; toc
    end
    
    % ------- Smooth strain field --------
    if DICpara.DispSmoothness>1e-6, USubpb2_Eulerian = funSmoothDispQuadtreeRBF(USubpb2_Eulerian,DICmesh,DICpara); end
    % ------- Don't smooth strain fields near the boundary --------
    for tempk=0:3, FSubpb2_Eulerian(4*DICmesh.markCoordHoleEdge-tempk) = FSubpb1_Eulerian(4*DICmesh.markCoordHoleEdge-tempk); end
    if DICpara.StrainSmoothness>1e-6, FSubpb2_Eulerian = funSmoothStrainQuadtreeRBF(0.1*FSubpb2_Eulerian+0.9*FSubpb1_Eulerian,DICmesh,DICpara); end
    % for tempk=0:1, USubpb2_Eulerian(2*markCoordHoleStrain-tempk) = USubpb1_Eulerian(2*markCoordHoleStrain-tempk); end
    % for tempk=0:3, FSubpb2_Eulerian(4*markCoordHoleStrain-tempk) = FSubpb1_Eulerian(4*markCoordHoleStrain-tempk); end
    for tempk=0:1, USubpb2_Eulerian(2*DICmesh.markCoordHoleEdge-tempk) = USubpb1_Eulerian(2*DICmesh.markCoordHoleEdge-tempk); end
    for tempk=0:3, FSubpb2_Eulerian(4*DICmesh.markCoordHoleEdge-tempk) = FSubpb1_Eulerian(4*DICmesh.markCoordHoleEdge-tempk); end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    % ------- Save data ------
    save(['Subpb2_step',num2str(ALSolveStep)],'USubpb2_Eulerian','FSubpb2_Eulerian');

    % ------ Plot ------
    USubpb2_EulerianWorld = USubpb2_Eulerian; USubpb2_EulerianWorld(2:2:end) = -USubpb2_Eulerian(2:2:end); FSubpb2_EulerianWorld = FSubpb2_Eulerian;  
    close all; Plotdisp_show(USubpb2_EulerianWorld,DICmesh.coordinatesFEM_EulerianWorld,DICmesh.elementsFEM(:,1:4),DICpara,'EdgeColor');
    Plotstrain_show(FSubpb2_EulerianWorld,DICmesh.coordinatesFEM_EulerianWorld,DICmesh.elementsFEM(:,1:4),DICpara,'EdgeColor');

    % ======= Update dual variables =======
    udual = FSubpb2_Eulerian - FSubpb1_Eulerian; vdual = USubpb2_Eulerian - USubpb1_Eulerian;
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
         
        %%%%% We use the fixed DIC subset window size %%%%%
        winsize_List = DICpara.winsize*ones(size(DICmesh.coordinatesFEM_Eulerian,1),2);
        DICpara.winsize_List = winsize_List;
         
        %%%%%%%%%%%%%%%%%%%%%%% Subproblem 1 %%%%%%%%%%%%%%%%%%%%%%%%%
        disp(['***** Start step',num2str(ALSolveStep),' Subproblem1 *****']);
        tic; [USubpb1_Eulerian,~,ALSub1Timetemp,ConvItPerEletemp,LocalICGNBadPtNumtemp] = Subpb1Quadtree(...
                                            USubpb2_Eulerian,FSubpb2_Eulerian,udual,vdual,DICmesh.coordinatesFEM_Eulerian,...
                                            Dg,gNormalized,fNormalized,mu,beta,HPar,ALSolveStep,DICpara,'GaussNewton',tol);
        FSubpb1_Eulerian = FSubpb2_Eulerian; toc 
        % for tempk=0:1, USubpb1_Eulerian(2*markCoordHoleStrain-tempk) = USubpb2_Eulerian(2*markCoordHoleStrain-tempk); end
        for tempk=0:1, USubpb1_Eulerian(2*DICmesh.markCoordHoleEdge-tempk) = USubpb2_Eulerian(2*DICmesh.markCoordHoleEdge-tempk); end
     
    
        ALSub1Time(ALSolveStep) = ALSub1Timetemp; ConvItPerEle(:,ALSolveStep) = ConvItPerEletemp; ALSub1BadPtNum(ALSolveStep) = LocalICGNBadPtNumtemp;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ------  Manually find some bad points from Local Subset ICGN step ------
        % disp('--- Start to manually remove bad points --- \n')
        % disp('    Comment codes here if you do not have bad local points \n')
        % %%%%% Comment START %%%%%
        %  [USubpb1_Eulerian,FSubpb1_Eulerian] = funRemoveOutliersQuadtree(DICmesh,DICpara,USubpb1_Eulerian,FSubpb1_Eulerian);
        %  disp('--- Remove bad points done ---')
        % %%%%% Comment END %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        save(['Subpb1_step',num2str(ALSolveStep)],'USubpb1_Eulerian','FSubpb1_Eulerian');
         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ============== Subproblem 2 ==============
        disp(['***** Start step',num2str(ALSolveStep),' Subproblem2 *****'])
        tic; [USubpb2_Eulerian] = Subpb2Quadtree(DICmesh,DICpara.GaussPtOrder,beta,mu,USubpb1_Eulerian,FSubpb1_Eulerian,udual,vdual,alpha,mean(DICpara.winstepsize),0);
		% [FSubpb2_Eulerian,~,~] = funGlobalNodalStrainQuadtree(DICmesh,USubpb2_Eulerian,DICpara.GaussPtOrder,0);
        [FSubpb2_Eulerian] = funGlobalNodalStrainRBF(DICmesh,DICpara,USubpb2_Eulerian);
        ALSub2Time(ALSolveStep) = toc; toc
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ------- Smooth strain field --------
        if DICpara.DispSmoothness>1e-6, USubpb2_Eulerian = funSmoothDispQuadtreeRBF(USubpb2_Eulerian,DICmesh,DICpara); end
        % ------- Don't change strain fields near the boundary --------
        for tempk=0:3, FSubpb2_Eulerian(4*DICmesh.markCoordHoleEdge-tempk) = FSubpb1_Eulerian(4*DICmesh.markCoordHoleEdge-tempk); end
        if DICpara.StrainSmoothness>1e-6, FSubpb2_Eulerian = funSmoothStrainQuadtreeRBF(0.1*FSubpb2_Eulerian+0.9*FSubpb1_Eulerian,DICmesh,DICpara); end
        % for tempk=0:1, USubpb2_Eulerian(2*markCoordHoleStrain-tempk) = USubpb1_Eulerian(2*markCoordHoleStrain-tempk); end
        % for tempk=0:3, FSubpb2_Eulerian(4*markCoordHoleStrain-tempk) = FSubpb1_Eulerian(4*markCoordHoleStrain-tempk); end
        for tempk=0:1, USubpb2_Eulerian(2*DICmesh.markCoordHoleEdge-tempk) = USubpb1_Eulerian(2*DICmesh.markCoordHoleEdge-tempk); end
        for tempk=0:3, FSubpb2_Eulerian(4*DICmesh.markCoordHoleEdge-tempk) = FSubpb1_Eulerian(4*DICmesh.markCoordHoleEdge-tempk); end
    
    
		save(['Subpb2_step',num2str(ALSolveStep)],'USubpb2_Eulerian','FSubpb2_Eulerian');
         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute norm of UpdateY
        USubpb2_Eulerian_Old = load(['Subpb2_step',num2str(ALSolveStep-1)],'USubpb2_Eulerian');
        USubpb2_Eulerian_New = load(['Subpb2_step',num2str(ALSolveStep)],'USubpb2_Eulerian');
        USubpb1_Eulerian_Old = load(['Subpb1_step',num2str(ALSolveStep-1)],'USubpb1_Eulerian');
        USubpb1_Eulerian_New = load(['Subpb1_step',num2str(ALSolveStep)],'USubpb1_Eulerian');
        if (mod(ImgSeqNum-2,DICpara.ImgSeqIncUnit) ~= 0 && (ImgSeqNum>2)) || (ImgSeqNum < DICpara.ImgSeqIncUnit)
            UpdateY = norm((USubpb2_Eulerian_Old.USubpb2_Eulerian - USubpb2_Eulerian_New.USubpb2_Eulerian), 2)/sqrt(size(USubpb2_Eulerian_Old.USubpb2_Eulerian,1));
            try
                UpdateY2 = norm((USubpb1_Eulerian_Old.USubpb1_Eulerian - USubpb1_Eulerian_New.USubpb1_Eulerian), 2)/sqrt(size(USubpb1_Eulerian_Old.USubpb1_Eulerian,1));
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
        udual = FSubpb2_Eulerian - FSubpb1_Eulerian; vdual = USubpb2_Eulerian - USubpb1_Eulerian; 
		 
        save(['uvdual_step',num2str(ALSolveStep)],'udual','vdual');
        try
        if UpdateY < tol2 || UpdateY2 < tol2
            break
        end
        catch
        end
         
    end
    fprintf('------------ Section 6 Done ------------ \n \n')
 
    
    
    %%%%% Convert into Lagrangian description %%%%%
    USubpb2 = -USubpb2_Eulerian;  
    
    temp_detF = FSubpb2_Eulerian(1:4:end) + FSubpb2_Eulerian(4:4:end) + ...
        FSubpb2_Eulerian(1:4:end).*FSubpb2_Eulerian(4:4:end) - ...
        FSubpb2_Eulerian(2:4:end).*FSubpb2_Eulerian(3:4:end) + 1;
    FSubpb2 = -FSubpb2_Eulerian;
    FSubpb2(1:4:end) = (FSubpb2_Eulerian(4:4:end)+1)./temp_detF - 1;
    FSubpb2(4:4:end) = (FSubpb2_Eulerian(1:4:end)+1)./temp_detF - 1;
    FSubpb2(2:4:end) = FSubpb2(2:4:end)./temp_detF;
    FSubpb2(3:4:end) = FSubpb2(3:4:end)./temp_detF;
    
    DICmesh.coordinatesFEM = DICmesh.coordinatesFEM_Eulerian - [USubpb2(1:2:end), USubpb2(2:2:end)];
 

 
    %%%%% Save data %%%%%
    ResultDisp{ImgSeqNum-1}.U = full(USubpb2);
    ResultDisp{ImgSeqNum-1}.U_Eulerian = full(USubpb2_Eulerian);
    ResultDisp{ImgSeqNum-1}.ALSub1BadPtNum = ALSub1BadPtNum;
    ResultDefGrad{ImgSeqNum-1}.F = full(FSubpb2);  
    ResultDefGrad{ImgSeqNum-1}.F_Eulerian = full(FSubpb2_Eulerian);  
    ResultFEMeshEachFrame{ImgSeqNum-1} = struct( 'coordinatesFEM_Eulerian',DICmesh.coordinatesFEM_Eulerian, ...
        'coordinatesFEM',DICmesh.coordinatesFEM,'elementsFEM',DICmesh.elementsFEM,'markCoordHoleEdge',DICmesh.markCoordHoleEdge );
    

%     % Save data
%     ResultDisp{ImgSeqNum-1}.U = full(USubpb1_Eulerian);
%     ResultDisp{ImgSeqNum-1}.ALSub1BadPtNum = ALSub1BadPtNum;
%     ResultDefGrad{ImgSeqNum-1}.F = full(FSubpb1_Eulerian); 
    

end    

%% ------ Plot ------
USubpb2_EulerianWorld = USubpb2_Eulerian; USubpb2_EulerianWorld(2:2:end) = -USubpb2_Eulerian(2:2:end); FSubpb2_EulerianWorld = FSubpb2_Eulerian; 
close all; Plotdisp_show(USubpb2_EulerianWorld,DICmesh.coordinatesFEM_EulerianWorld,DICmesh.elementsFEM(:,1:4),DICpara,'EdgeColor');
Plotstrain_show(FSubpb2_EulerianWorld,DICmesh.coordinatesFEM_EulerianWorld,DICmesh.elementsFEM(:,1:4),DICpara,'EdgeColor');
 
% winsize_xy = [winsize_x_List, winsize_y_List]'; winsize_xy = winsize_xy(:);
% Plotdisp_show(winsize_xy,DICmesh.coordinatesFEM_EulerianWorld,DICmesh.elementsFEM(:,1:4),DICpara,'EdgeColor');

% figure, plot3(DICmesh.coordinatesFEM_Eulerian(:,1),DICmesh.coordinatesFEM_Eulerian(:,2), USubpb2_Eulerian(2:2:end),'.');


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
    USubpb2_Eulerian = load(['Subpb2_step',num2str(ALSolveStep )],'USubpb2_Eulerian');
    USubpb1_Eulerian = load(['Subpb1_step',num2str(ALSolveStep )],'USubpb1_Eulerian');
    UpdateY = norm((USubpb2_Eulerian.USubpb2_Eulerian - USubpb1_Eulerian.USubpb1_Eulerian), 2)/sqrt(length(USubpb2_Eulerian.USubpb2_Eulerian));
    disp(num2str(UpdateY));
end
disp('==== Fhat^(k) - F^(k) ====');
for ALSolveStep = 1:ALSolveStep1
    FSubpb1_Eulerian = load(['Subpb1_step',num2str(ALSolveStep )],'FSubpb1_Eulerian');
    FSubpb2_Eulerian = load(['Subpb2_step',num2str(ALSolveStep )],'FSubpb2_Eulerian');
    UpdateF = norm((FSubpb1_Eulerian.FSubpb1_Eulerian - FSubpb2_Eulerian.FSubpb2_Eulerian), 2)/sqrt(length(FSubpb1_Eulerian.FSubpb1_Eulerian));
    disp(num2str(UpdateF));
end
disp('==== uhat^(k) - uhat^(k-1) ====');
for ALSolveStep = 2:ALSolveStep1
    USubpb2_Eulerian_Old = load(['Subpb2_step',num2str(ALSolveStep-1)],'USubpb2_Eulerian');
    USubpb2_Eulerian_New = load(['Subpb2_step',num2str(ALSolveStep)],'USubpb2_Eulerian');
    UpdateY = norm((USubpb2_Eulerian_Old.USubpb2_Eulerian - USubpb2_Eulerian_New.USubpb2_Eulerian), 2)/sqrt(length(USubpb2_Eulerian.USubpb2_Eulerian));
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
clear coordinatesFEM_EulerianQuadtree elementsFEMQuadtree 


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

%% ------ Start main part ------
for ImgSeqNum = 2 % [2 : length(ImgNormalized)] 
    
    close all; 
    disp(['Current image frame #: ', num2str(ImgSeqNum),'/',num2str(length(ImgNormalized))]);
    DICpara.ImgRefMask = fNormalizedMask;
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    USubpb2 = ResultDisp{ImgSeqNum-1}.U;
    USubpb2_Eulerian = -USubpb2;
    FSubpb2 = ResultDefGrad{ImgSeqNum-1}.F;
    try % Just in case you haven't saved "F_Eulerian" before
        FSubpb2_Eulerian = ResultDefGrad{ImgSeqNum-1}.F_Eulerian;
    catch
        temp_detF = FSubpb2(1:4:end) + FSubpb2(4:4:end) + ...
        FSubpb2(1:4:end).*FSubpb2(4:4:end) - ...
        FSubpb2(2:4:end).*FSubpb2(3:4:end) + 1;
        FSubpb2_Eulerian = -FSubpb2;
        FSubpb2_Eulerian(1:4:end) = (FSubpb2(4:4:end)+1)./temp_detF - 1;
        FSubpb2_Eulerian(4:4:end) = (FSubpb2(1:4:end)+1)./temp_detF - 1;
        FSubpb2_Eulerian(2:4:end) = FSubpb2_Eulerian(2:4:end)./temp_detF;
        FSubpb2_Eulerian(3:4:end) = FSubpb2_Eulerian(3:4:end)./temp_detF;
    end
    coordinatesFEM = ResultFEMeshEachFrame{ImgSeqNum-1}.coordinatesFEM;
    coordinatesFEM_Eulerian = ResultFEMeshEachFrame{ImgSeqNum-1}.coordinatesFEM_Eulerian;
    elementsFEM = ResultFEMeshEachFrame{ImgSeqNum-1}.elementsFEM;
    try markCoordHoleEdge = ResultFEMeshEachFrame{ImgSeqNum-1}.markCoordHoleEdge; catch; end
    DICmesh.coordinatesFEM = coordinatesFEM;
    DICmesh.elementsFEM = elementsFEM;
    coordinatesFEMWorld = DICpara.um2px*[coordinatesFEM(:,1),size(ImgNormalized{1},2)+1-coordinatesFEM(:,2)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    % ------ Plotting and Compute Strain-------
    if size(USubpb2,1) == 1
        ULocal = USubpb2_New.USubpb2; FLocal = FSubpb2.FSubpb2; 
    else
        ULocal = USubpb2; FLocal = FSubpb2; ULocal_Eulerian = USubpb2_Eulerian; 
    end
     
    % ------ Smooth displacements ------
    %prompt = 'Do you want to smooth displacement? (0-yes; 1-no)';
    %DoYouWantToSmoothOnceMore = input(prompt); 
    SmoothTimes = 0;
    try
        while DICpara.DoYouWantToSmoothOnceMore == 0 && SmoothTimes < 1
            ULocal = funSmoothDispQuadtreeRBF(ULocal,DICmesh,DICpara);
            %close all; Plotuv(ULocal,x0,y0); %DICpara.DoYouWantToSmoothOnceMore = input(prompt);
            SmoothTimes = SmoothTimes + 1;
        end
    catch
    end
    
     
    % ----- Compute strain field ------
    %%%%% Load deformed image %%%%%%%
    gNormalized = ImgNormalized{ ImgSeqNum } ; % Load current deformed image frame 
    gNormalizedMask = double( ImgMask{ImgSeqNum} );
    gNormalized = gNormalized .* gNormalizedMask ;
    Dg = funImgGradient(gNormalized,gNormalized,gNormalizedMask); % Finite difference to compute image grayscale gradients;
    
    %%%%% Compute strain field %%%%%
    ComputeStrainQuadtree;  
    
    %%%%% Scale to the physical world %%%%%
    UWorld = DICpara.um2px*ULocal; UWorld(2:2:end) = -UWorld(2:2:end);   % close all; Plotuv(UWorld,x0,y0World);
      
    
    % ------ Plot disp and strain ------
    if DICpara.OrigDICImgTransparency == 1
        Plotdisp_show(UWorld,coordinatesFEMWorld,DICmesh.elementsFEM(:,1:4),DICpara,'NoEdgeColor');
        [strain_exx,strain_exy,strain_eyy,strain_principal_max,strain_principal_min,strain_maxshear,strain_vonMises] = ...
                   Plotstrain0Quadtree(FStraintemp,coordinatesFEMWorld,elementsFEM(:,1:4),DICpara);
    
    else % Plot over raw DIC images
        if DICpara.Image2PlotResults == 0 % Plot over the first image; "file_name{1,1}" corresponds to the first image	
             
            % PlotdispQuadtree(UWorld,coordinatesFEMWorld,elementsFEM(:,1:4),...
            %     file_name{1,1}, DICpara);
            % [strain_exx,strain_exy,strain_eyy,strain_principal_max,strain_principal_min, ...
            %     strain_maxshear,strain_vonMises] = PlotstrainQuadtree(UWorld,FStraintemp, ...
            %     coordinatesFEMWorld,elementsFEM(:,1:4),file_name{1,1},DICpara);
            PlotdispQuadtreeMasks(UWorld,coordinatesFEMWorld,elementsFEM(:,1:4),...
                file_name{1,1}, ImgMask{1},DICpara);
             
            [strain_exx,strain_exy,strain_eyy,strain_principal_max,strain_principal_min, ...
               strain_maxshear,strain_vonMises] = PlotstrainQuadtreeMasks(UWorld,FStraintemp, ...
               coordinatesFEMWorld,elementsFEM(:,1:4),file_name{1,1}, ImgMask{1},DICpara);
           
           
        else % Plot over second or next deformed images
             
             PlotdispQuadtreeMasks(UWorld,coordinatesFEMWorld,elementsFEM(:,1:4),...
                file_name{1, ImgSeqNum}, ImgMask{ ImgSeqNum },DICpara);
            
            [strain_exx,strain_exy,strain_eyy,strain_principal_max,strain_principal_min, ...
               strain_maxshear,strain_vonMises] = PlotstrainQuadtreeMasks(UWorld,FStraintemp, ...
               coordinatesFEMWorld,elementsFEM(:,1:4),file_name{1, ImgSeqNum}, ...
               ImgMask{ ImgSeqNum },DICpara);
             
        end
    end
 
    % ----- Save strain results ------
    ResultStrain{ImgSeqNum-1} = struct('strainxCoord',coordinatesFEMWorld(:,1),'strainyCoord',coordinatesFEMWorld(:,2), ...
            'dispu',UWorld(1:2:end),'dispv',UWorld(2:2:end), ...
            'dudx',FStraintemp(1:4:end),'dvdx',FStraintemp(2:4:end),'dudy',FStraintemp(3:4:end),'dvdy',FStraintemp(4:4:end), ...
            'strain_exx',strain_exx,'strain_exy',strain_exy,'strain_eyy',strain_eyy, ...
            'strain_principal_max',strain_principal_max,'strain_principal_min',strain_principal_min, ...
            'strain_maxshear',strain_maxshear,'strain_vonMises',strain_vonMises);
      
    % ------ Save figures for tracked displacement and strain fields ------
    % SaveFigFilesDispAndStrainQuadtree;
    % pause;
    
    
end
% ------ END of for-loop {ImgSeqNum = 2:length(ImgNormalized)} ------
fprintf('------------ Section 8 Done ------------ \n \n')


% ------ Save data again including solved strain fields ------
% results_name = ['results_',imgname,'_ws',num2str(DICpara.winsize),'_st',num2str(DICpara.winstepsize),'.mat'];
% save(results_name, 'file_name','DICpara','DICmesh','ResultDisp','ResultDefGrad','ResultFEMesh','ResultFEMeshEachFrame',...
%                    'ALSub1Time','ALSub2Time','ALSolveStep','ResultStrain');
% 


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
v = VideoWriter('video_mesh.mp4');
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

