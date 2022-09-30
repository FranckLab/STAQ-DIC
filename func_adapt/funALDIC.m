
function [U0,ULoc,FLoc,USubpb2,FSubpb2,udual,vdual,ConvItPerEle,ALSolveStep,coordinatesFEM,elementsFEM,dirichlet,neumann,mu,beta,...
    ALSub1Time,ALSub2Time] = funALDIC(ImgNormalized,file_name,ImgSeqNum,gridxyROIRange,...
    winsize,winstepsize,ClusterNo,Subpb2FDOrFEM)


    %% Section 3
    fprintf('\n'); fprintf('------------ Section 3 Start ------------ \n')
    % % Do first deformed image at first.
    gridxROIRange = gridxyROIRange.gridx; gridyROIRange = gridxyROIRange.gridy; 
    fNormalized = ImgNormalized{1}; gNormalized = ImgNormalized{ImgSeqNum};
    
    if ImgSeqNum == 2 || StillFFTSearch == 1
        % ====== Integer Search ======
        [SizeOfFFTSearchRegion,x0,y0,u,v,cc]= IntegerSearch(gridxROIRange,gridyROIRange,fNormalized,gNormalized,winsize,winstepsize,file_name);
        % ====== FEM mesh set up ======
        [coordinatesFEM,elementsFEM,coordinates,elements,dirichlet,neumann,x0,y0] = MeshSetUp(x0,y0,winstepsize);
        % ====== Remove outliers ======
        % Already included in the above IntegerSearch.m function
        % qDICOrNot = 0; Thr0 = 100; [u,v,cc] = funRemoveOutliers(u,v,cc,qDICOrNot,Thr0);
        % ====== Initial Value ======
        U0 = Init(u,v,cc.max,x0,y0,0); PlotuvInit; 
    else
        U0 = ResultDisp{ImgSeqNum-2}.U;
    end
    
    %Plotdisp_show(U0,[coordinatesFEM(:,1),size(fNormalized,2)+1-coordinatesFEM(:,2)],elementsFEM); % Plot initial values
    % ====== Spline interpolation images ======
    %[imgfNormalizedbc,imggNormalizedbc,imgSize,DfAxis] = funImgGradient(fNormalized,gNormalized);
    Df = funImgGradient(fNormalized,gNormalized); % % using finite difference;
    % ====== Compute f(X)-g(x+u) ======
    % PlotImgDiff(x0,y0,u,v,fNormalized,gNormalized);
    fprintf('------------ Section 3 Done ------------ \n \n')


    %% Section 4
    fprintf('------------ Section 4 Start ------------ \n')
    % ====== ALStep 1 Subproblem1: Local Subset DIC ======
    mu = 0; beta = 0; tol = 1e-6; ALSolveStep = 1; ALSub1Time = zeros(6,1); ALSub2Time = zeros(6,1); ConvItPerEle = zeros(size(coordinatesFEM,1),6);
    disp(['***** Start step',num2str(ALSolveStep),' Subproblem1 *****'])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ------ Assign parpool cluster No ------
    % prompt = 'How many parallel pools to open? (Put in 1 if no parallel computing): ';  ClusterNo = input(prompt);
    % ------ Start Local DIC IC-GN iteration ------
    [USubpb1,FSubpb1,HtempPar,ALSub1Timetemp,ConvItPerEletemp] = LocalICGN(U0,coordinatesFEM,...
        Df,fNormalized,gNormalized,winsize,winstepsize,tol,'GaussNewton',ClusterNo);
    ALSub1Time(ALSolveStep) = ALSub1Timetemp; ConvItPerEle(:,ALSolveStep) = ConvItPerEletemp; clear ALSub1Timetemp ConvItPerEletemp; toc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ------ Save local subset DIC results ------
    ULoc = USubpb1; FLoc = FSubpb1;
    % ------ Manually find some bad points from Local Subset ICGN step ------
    disp('--- Start to manually remove bad points ---')
    USubpb1World = USubpb1; USubpb1World(2:2:end) = -USubpb1(2:2:end); Plotuv(USubpb1World,x0,y0World);
    u = reshape(USubpb1(1:2:end),size(x0,1),size(x0,2)); v = reshape(USubpb1(2:2:end),size(x0,1),size(x0,2));
    qDICOrNot = 0; Thr0 = 100; [u,v,~,Local_BadptRow,Local_BadptCol,RemoveOutliersList] = funRemoveOutliers(u',v',[],qDICOrNot,Thr0); u=u';v=v';
    USubpb1(1:2:end) = reshape(u,size(elements,1),1); USubpb1(2:2:end) = reshape(v,size(elements,1),1);
    f11 = reshape(FSubpb1(1:4:end),size(x0,1),size(x0,2)); f21 = reshape(FSubpb1(2:4:end),size(x0,1),size(x0,2)); 
    f12 = reshape(FSubpb1(3:4:end),size(x0,1),size(x0,2)); f22 = reshape(FSubpb1(4:4:end),size(x0,1),size(x0,2)); 
    f11=f11'; f21=f21'; f12=f12'; f22=f22';
    f11(RemoveOutliersList) = NaN; f21(RemoveOutliersList) = NaN; f12(RemoveOutliersList) = NaN; f22(RemoveOutliersList) = NaN;
    f11 = inpaint_nans(f11,4); f21 = inpaint_nans(f21,4); f12 = inpaint_nans(f12,4); f22 = inpaint_nans(f22,4);
    f11=f11'; f21=f21'; f12=f12'; f22=f22';
    FSubpb1(1:4:end) = reshape(f11,size(elements,1),1); FSubpb1(2:4:end) = reshape(f21,size(elements,1),1);
    FSubpb1(3:4:end) = reshape(f12,size(elements,1),1); FSubpb1(4:4:end) = reshape(f22,size(elements,1),1);
    disp('--- Remove bad points done ---')
    % ------ Plot ------
    USubpb1World = USubpb1; USubpb1World(2:2:end) = -USubpb1(2:2:end); 
    FSubpb1World = FSubpb1; % FSubpb1World(2:2:end) = -FSubpb1(2:2:end);
    close all; Plotuv(USubpb1World,x0,y0World); Plotdisp_show(USubpb1World,coordinatesFEMWorld,elementsFEM);
    Plotstrain_show(FSubpb1World,coordinatesFEMWorld,elementsFEM);
    save(['Subpb1_step',num2str(ALSolveStep)],'USubpb1','FSubpb1');
    fprintf('------------ Section 4 Done ------------ \n \n')


    %% Section 5
    fprintf('------------ Section 5 Start ------------ \n')
    % ====== Choose solver method ======
    try
        if Subpb2FDOrFEM == 0  %Algorithm to solve Subproblem 2. 1): Finite difference (DF); 2): Finite element (FE);
            fprintf('Method to solve ALDIC Subproblem 2:    \n')
            fprintf('1: Finite difference   \n')
            fprintf('2: Finite element method  \n')
            prompt = 'Input here: ';
            Subpb2FDOrFEM = input(prompt);
        end
    catch
        Subpb2FDOrFEM = 1;
    end
    
    % ======= ALStep 1 Subproblem 2: Global constraint =======
    % ------ Smooth displacements for better F ------
    DispFilterSize=0; DispFilterStd = 0; StrainFilterSize = 0; StrainFilterStd =0; LevelNo = 1;
    FSubpb1 = funSmoothStrain(FSubpb1,coordinatesFEM,elementsFEM,winstepsize,StrainFilterSize,StrainFilterStd);
    % % ------ Smooth displacements ------
    % prompt = 'Do you want to smooth displacement? (0-yes; 1-no)';
    % DoYouWantToSmoothOnceMore = input(prompt); DispFilterSize=0; DispFilterStd=0;
    % if DoYouWantToSmoothOnceMore == 0,USubpb1 = funSmoothDisp(USubpb1,coordinatesFEM,elementsFEM,winstepsize,DispFilterSize,DispFilterStd);end

	% ====== Define penalty parameter ======
    mu = 1e-3; udual = 0*FSubpb1; vdual = 0*USubpb1; 
    betaList = [1e-3,sqrt(1e-5),1e-2,sqrt(1e-3),1e-1,sqrt(1e-1)]*mean(winstepsize).^2.*mu; 
    Err1 = zeros(length(betaList),1); Err2 = Err1;
    
    disp(['***** Start step',num2str(ALSolveStep),' Subproblem2 *****'])
    % ====== Decide to use FD or FE methods to solve Subpb2 step ======
	if Subpb2FDOrFEM == 1 % Using FD method
		% ====== Build sparse finite difference operator ======
		disp('Assemble finite difference operator D');
		M = size(x0,1); N = size(x0,2);
		tic; Rad = 1; D = funDerivativeOp((M-2*Rad),(N-2*Rad),winstepsize); % D = sparse(4*(M-2*Rad)*(N-2*Rad), 2*(M-2*Rad)*(N-2*Rad));
		D2 = funDerivativeOp(M,N,winstepsize); toc
		disp('Finish assembling finite difference operator D');
		% ===== Solver using finite difference approximation ======
		tic; a = FSubpb1-udual; b = USubpb1-vdual; 
		Rad = 1; FDNeumannBC; % Find coordinatesFEM that belong to (x(Rad+1:M-Rad,Rad+1:N-Rad),y(Rad+1:M-Rad,Rad+1:N-Rad))
		atemp = a(temp3); btemp = b(temp4);
        for tempk = 1:length(betaList)
            beta = betaList(tempk);
            tempAMatrixSub2 = (beta*(D')*D) + mu*speye(2*(M-2*Rad)*(N-2*Rad));
            USubpb2temp = (tempAMatrixSub2) \ (beta*D'*atemp + mu*btemp ) ;
            USubpb2 = USubpb1; USubpb2(temp4) = USubpb2temp; 
            FSubpb2 = D2*USubpb2; 
            
            Err1(tempk) = norm(USubpb1-USubpb2,2);
            Err2(tempk) = norm(FSubpb1-FSubpb2,2);
        end
        ErrSum = Err1+Err2*mean(winstepsize)^2; [~,indexOfbeta] = min(ErrSum);
        try
            [fitobj] = fit(log10(betaList(indexOfbeta-1:1:indexOfbeta+1))',ErrSum(indexOfbeta-1:1:indexOfbeta+1),'poly2');
            p = coeffvalues(fitobj); beta = 10^(-p(2)/2/p(1));
        catch
            beta = betaList(indexOfbeta);
        end
    %%JY!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %if Subpb2FDOrFEM == 1 % Using FD method
    %%JY!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Using optimal beta to solve again
        tempAMatrixSub2 = (beta*(D')*D) + mu*speye(2*(M-2*Rad)*(N-2*Rad));
        USubpb2temp = (tempAMatrixSub2) \ (beta*D'*atemp + mu*btemp) ;
        USubpb2 = USubpb1; USubpb2(temp4) = USubpb2temp;
		%%%%%%%%%%%%%% End of using finite difference approximation %%%%%%%%%%%%%%
	else % Subpb2FDOrFEM == 2: Using FE method
         M = size(x0,1); N = size(x0,2); GaussPtOrder = 2; alpha = 0;
        close all;  
        % ====== Solver using finite element method ======
        for tempk = 1:length(betaList)
            beta = betaList(tempk);
            GaussPtOrder = 2; alpha = 0; [USubpb2] = Subpb2(coordinatesFEM,elementsFEM,M,N,neumann,beta,mu,USubpb1,FSubpb1,udual,vdual,alpha,GaussPtOrder,ClusterNo);
            [FSubpb2] = funGlobal_NodalStrainAvg(coordinatesFEM,elementsFEM,USubpb2,GaussPtOrder);
            % Plotstrain_show(FSubpb2,coordinatesFEM,elementsFEM);
            % Plotstrain_show(FSubpb1,coordinatesFEM,elementsFEM);
            Err1(tempk) = norm(USubpb1-USubpb2,2);
            Err2(tempk) = norm(FSubpb1-FSubpb2,2);
        end
        Err1Norm = (Err1-mean(Err1))/std(Err1); figure, plot(Err1Norm);
        Err2Norm = (Err2-mean(Err2))/std(Err2); figure, plot(Err2Norm);
        ErrSum = Err1Norm+Err2Norm; figure,plot(ErrSum); [~,indexOfbeta] = min(ErrSum);
        try
            [fitobj] = fit(log10(betaList(indexOfbeta-1:1:indexOfbeta+1))',ErrSum(indexOfbeta-1:1:indexOfbeta+1),'poly2');
            p = coeffvalues(fitobj); beta = 10^(-p(2)/2/p(1));
        catch
            beta = betaList(indexOfbeta);
        end
        % Using optimal beta to solve again
        [USubpb2] = Subpb2(coordinatesFEM,elementsFEM,M,N,neumann,beta,mu,USubpb1,FSubpb1,udual,vdual,alpha,GaussPtOrder,ClusterNo);
    end
	ALSub2Time(ALSolveStep) = toc; toc
	
    % ------- Before computing strain, we smooth the displacement field -------
    %USubpb2 = funSmoothDisp(USubpb2,coordinatesFEM,elementsFEM,x0,y0,winstepsize,DispFilterSize,DispFilterStd);
    % ------- Compute strain field --------
    if Subpb2FDOrFEM == 1 %FD 
		FSubpb2 = D2*USubpb2; % D2 = funDerivativeOp(M,N,winstepsize);
	else %FEM
		[FSubpb2] = funGlobal_NodalStrainAvg(coordinatesFEM,elementsFEM,USubpb2,GaussPtOrder);
    end
    
    % ------- Smooth strain field --------
    % FSubpb2 = funSmoothStrain(FSubpb2,coordinatesFEM,elementsFEM,winstepsize,StrainFilterSize,StrainFilterStd);

    % ------- Save data ------
    save(['Subpb2_step',num2str(ALSolveStep)],'USubpb2','FSubpb2');

    % ------ Plot ------
    USubpb2World = full(USubpb2); USubpb2World(2:2:end) = (-USubpb2(2:2:end)); 
    FSubpb2World = full(FSubpb2); % FSubpb1World(2:2:end) = -FSubpb1(2:2:end);
    close all; %Plotuv(USubpb2World,x0,y0World); 
    Plotdisp_show(USubpb2World,coordinatesFEMWorld,elementsFEM);
    Plotstrain_show(FSubpb2World,coordinatesFEMWorld,elementsFEM);

    % ======= Update dual variables =======
    if Subpb2FDOrFEM == 1 %FD
		udualtemp1 = (FSubpb2 - FSubpb1); udualtemp2 = udualtemp1(temp3);
		vdualtemp1 = (USubpb2 - USubpb1); vdualtemp2 = vdualtemp1(temp4);
		udual = zeros(4*M*N,1); vdual = zeros(2*M*N,1);
		udual(temp3) = udualtemp2; vdual(temp4) = vdualtemp2;
	else  % FEM or other methods
		udual = FSubpb2 - FSubpb1; vdual = USubpb2 - USubpb1;
	end
    save(['uvdual_step',num2str(ALSolveStep)],'udual','vdual');
    fprintf('------------ Section 5 Done ------------ \n \n')
 
    %Plotstrain_show(FSubpb2,coordinatesFEM,elementsFEM);
   
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Section 6
    fprintf('------------ Section 6 Start ------------ \n')
    % ==================== Global AL Loop ==========================
    ALSolveStep = 1; tol2 = 1e-4; UpdateY = 1e4; CrackOrNot = 0; CrackPath1 = [0,0]; CrackPath2 = [0,0]; CrackTip = [0,0]; 
    HPar = cell(21,1); for tempj = 1:21, HPar{tempj} = HtempPar(:,tempj); end

    while (ALSolveStep < 3)
        ALSolveStep = ALSolveStep + 1;  % Update using the last step
        %%%%%%%%%%%%%%%%%%%%%%% Subproblem 1 %%%%%%%%%%%%%%%%%%%%%%%%%
        disp(['***** Start step',num2str(ALSolveStep),' Subproblem1 *****'])
        tic;[USubpb1, HPar, ALSub1Timetemp, ConvItPerEletemp] = Subpb1(USubpb2,FSubpb2,udual,vdual,coordinatesFEM,...
        Df,fNormalized,gNormalized,winsize,mu,beta,HPar,ALSolveStep,tol,'GaussNewton',ClusterNo);
        USubpb1 = full(USubpb1); FSubpb1 = full(FSubpb2); toc
        ALSub1Time(ALSolveStep) = ALSub1Timetemp; ConvItPerEle(:,ALSolveStep) = ConvItPerEletemp; clear ALSub1Timetemp ConvItPerEletemp;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ------  Manually find some bad points from Local Subset ICGN step ------
        disp('--- Start to manually remove bad points ---')
        close all; Plotuv(USubpb1,x0,y0);
        u = reshape(USubpb1(1:2:end),M,N); v = reshape(USubpb1(2:2:end),M,N);
        %[u,v,~,Subpb1_BadptRow,Subpb1_BadptCol] = funRemoveOutliers(u,v,[],0,Local_BadptRow,Local_BadptCol); 
        [u,v,cc] = funRemoveOutliers(u,v,cc,qDICOrNot,Thr0,Local_BadptRow,Local_BadptCol);
        disp('--- Remove bad points done ---')
        USubpb1(1:2:end) = reshape(u,size(elements,1),1); USubpb1(2:2:end) = reshape(v,size(elements,1),1);
        close all; Plotuv(USubpb1,x0,y0); Plotdisp_show(USubpb1,coordinatesFEM,elementsFEM);
        save(['Subpb1_step',num2str(ALSolveStep)],'USubpb1','FSubpb1');
        % %USubpb1 = funSmoothDisp(USubpb1,coordinatesFEM,elementsFEM,x0,y0,winstepsize,DispFilterStd,DispFilterSize,LevelNo);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ============== Subproblem 2 ==============
        disp(['***** Start step',num2str(ALSolveStep),' Subproblem2 *****'])
        if Subpb2FDOrFEM == 1 %FD
			% ------- using finite difference approximation --------
			tic; a = FSubpb1-udual; b = USubpb1-vdual; atemp = a(temp3); btemp = b(temp4);
			USubpb2temp = (tempAMatrixSub2) \ (beta*D'*atemp + mu*btemp ) ;
            USubpb2 = USubpb1; USubpb2(temp4) = USubpb2temp; %toc
			% ------- End of using finite difference approximation --------
		else % FEM 
			tic; [USubpb2] = Subpb2(coordinatesFEM,elementsFEM,M,N,neumann,beta,mu,USubpb1,FSubpb1,udual,vdual,alpha,GaussPtOrder,ClusterNo);
        end
		ALSub2Time(ALSolveStep) = toc; toc
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % ------- Before computing strain, we smooth the displacement field -------
        %USubpb2 = funSmoothDisp(USubpb2,coordinatesFEM,elementsFEM,x0,y0,winstepsize,DispFilterSize,DispFilterStd);
        % ------- Compute strain field --------
        if Subpb2FDOrFEM == 1 %FD
			FSubpb2 = D2*USubpb2; % D2 = funDerivativeOp(M,N,winstepsize);
		else %FEM
			GaussPtOrder = 2; [FSubpb2] = funGlobal_NodalStrainAvg(coordinatesFEM,elementsFEM,USubpb2,GaussPtOrder);
		end
		
		% ------- Smooth strain field --------
        %FSubpb2 = funSmoothStrain(FSubpb2,coordinatesFEM,elementsFEM,winstepsize,StrainFilterSize,StrainFilterStd);
        %Plotstrain_show(FSubpb2,coordinatesFEM,elementsFEM);
        save(['Subpb2_step',num2str(ALSolveStep)],'USubpb2','FSubpb2');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute norm of UpdateY
        USubpb2_Old = load(['Subpb2_step',num2str(ALSolveStep-1)],'USubpb2');
        USubpb2_New = load(['Subpb2_step',num2str(ALSolveStep)],'USubpb2');
        UpdateY = norm((USubpb2_Old.USubpb2 - USubpb2_New.USubpb2), 2)/sqrt(size(USubpb2_Old.USubpb2,1));

        USubpb1_Old = load(['Subpb1_step',num2str(ALSolveStep-1)],'USubpb1');
        USubpb1_New = load(['Subpb1_step',num2str(ALSolveStep)],'USubpb1');
        UpdateY2 = norm((USubpb1_Old.USubpb1 - USubpb1_New.USubpb1), 2)/sqrt(size(USubpb1_Old.USubpb1,1));

        disp(['Update local step  = ',num2str(UpdateY2)]);
        disp(['Update global step = ',num2str(UpdateY)]);
        fprintf('*********************************** \n \n');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Update dual variables------------------------------
        if Subpb2FDOrFEM == 1 %FD
			udualtemp1 =  (FSubpb2 - FSubpb1); udualtemp2 = udualtemp1(temp3);
			vdualtemp1 =  (USubpb2 - USubpb1); vdualtemp2 = vdualtemp1(temp4);
			udual(temp3) = udual(temp3)+udualtemp2; 
			vdual(temp4) = vdual(temp4)+vdualtemp2;
		else %FEM
			udual = FSubpb2 - FSubpb1; vdual = USubpb2 - USubpb1; 
		end

        save(['uvdual_step',num2str(ALSolveStep)],'udual','vdual');

        if UpdateY < tol2 || UpdateY2 < tol2
            break
        end

    end
    fprintf('------------ Section 6 Done ------------ \n \n')

  
    
    % ------ Delete temp files ------
    for tempi = 1:ALSolveStep
        file_name_Subpb1 = ['Subpb1_step',num2str(tempi),'.mat'];
        file_name_Subpb2 = ['Subpb2_step',num2str(tempi),'.mat'];
        file_name_dual = ['uvdual_step',num2str(tempi),'.mat'];
        delete(file_name_Subpb1); delete(file_name_Subpb2); delete(file_name_dual);
    end
    