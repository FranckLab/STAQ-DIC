%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that do adaptive global DIC with Kuhn triangulation mesh
% in the oder of: (Input: (LevelNo-1)Solve) -> Estimate -> Mark -> Refine -> Solve(Output)
% 
% In this code, variables at (LevelNo-1) are all called "Level1", while
% variables at (LevelNo) are called "Level2".
%
% Author: Jin Yang, Email: jyang526@wisc.edu
% Date: 2019 June.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ULevel2,FLevel2,alphavar,normOfWLevel2,...
    coordinatesFEMLevel2,elementsFEMLevel2,eleGenerationLevel2,dirichletLevel2,neumannLevel2,...
    thetaDorfler,TimeICGNLevel2,EstimateTime,MarkTime,RefineTime,...
    rhoLevel1Vector,rhoLevel1Vectortemp1,rhoLevel1Vectortemp2,EnrHAndTipEleIndex,EnrTipEleIndex] ...
    = funGBMSIterTri(ULevel1,FLevel1,coordinatesFEM,elementsFEM,winstepsize,ClusterNo,...
    eleGeneration,dirichlet,neumann,LevelNo,fNormalized,gNormalized,winsize,...
    CrackOrNot,CrackPath1,CrackPath2,CrackTip,CrackTipOrNot,EnrHAndTipEleIndex,EnrTipEleIndex,...
    DispFilterSize,DispFilterStd,StrainFilterSize,StrainFilterStd,tol,alphavarLevel1)

% ULevel1 = UIter; FLevel1 = FIter; udualLevel1 = udualIter; vdualLevel1 = vdualIter; USubpb2tempLevel1 = USubpb2tempIter;
% elementsFEM = elementsFEMIter;  coordinatesFEM = coordinatesFEMIter;
% eleGeneration = eleGenerationIter; dirichlet = dirichletIter; neumann = neumannIter;
% fNormalized = ImgNormalized{1}; gNormalized = ImgNormalized{2};
% betamuLevel1 = betamuIter;
alphavar = alphavarLevel1; GaussPtOrder = 1;
imgPyramidUnit = 1; Df = funImgGradient(fNormalized,gNormalized); 

%% ============== Estimate: Level1 ==============
% ------ Compute rho to find where to refine ------
rhoLevel1Vector = zeros(size(elementsFEM,1),1); rhoLevel1Vectortemp1 = rhoLevel1Vector; rhoLevel1Vectortemp2 = rhoLevel1Vector;
imgPyramidUnit = 1; hbar = parfor_progressbar(size(elementsFEM,1),'Please wait for a posterior error estimate!'); tic

% ------ Compute CST strain Nabla_U ------
[~,CSTStrain,~] = Global_CSTStrain(coordinatesFEM,elementsFEM,ULevel1);

parfor j = 1:size(elementsFEM,1)
    
    % Find element node indices
    temp = [2*elementsFEM(j,1)-1, 2*elementsFEM(j,1), 2*elementsFEM(j,2)-1, 2*elementsFEM(j,2), 2*elementsFEM(j,3)-1, 2*elementsFEM(j,3)];
    % Find eleCase of the triangle element: 8 different orientations
    [longestEdgeMidPtCoord,longestEdgeNo,longestEdgeLength,eleCase] = funFindLongestEdge(coordinatesFEM(elementsFEM(j,1:3),:));
    % Find neighboring elements for jump residuals
    if (eleCase ~= 0)
        [eleNeighbor2,eleNeighbor3,eleNeighbor4] = funFindEleNeighbors(elementsFEM,j,longestEdgeNo);
    else
        disp('There is something wrong with the findLongestEdge!')
    end
    
    CSTNeighborAtBoundaryCheck = zeros(4,1); CSTStraintemp = zeros(4,4); CSTStraintemp(1,1:4) = CSTStrain(j,1:4);
    if eleNeighbor2 ~= 0, CSTStraintemp(2,1:4) = CSTStrain(eleNeighbor2,1:4); CSTNeighborAtBoundaryCheck(2) = 1; end
    if eleNeighbor3 ~= 0, CSTStraintemp(3,1:4) = CSTStrain(eleNeighbor3,1:4); CSTNeighborAtBoundaryCheck(3) = 1; end
    if eleNeighbor4 ~= 0, CSTStraintemp(4,1:4) = CSTStrain(eleNeighbor4,1:4); CSTNeighborAtBoundaryCheck(4) = 1; end
    
    % ---- Compute error estimator ----
    rhoLevel1Vectortemp1(j) = (longestEdgeLength/winstepsize)^2 * aPostErrEstIntTri(coordinatesFEM(elementsFEM(j,:),:),fNormalized,gNormalized,Df,ULevel1(temp)) ;
    rhoLevel1Vectortemp2(j) = (longestEdgeLength/winstepsize) * alphavar * aPostErrEstJumpTri(eleCase,longestEdgeLength,CSTStraintemp,CSTNeighborAtBoundaryCheck);
    rhoLevel1Vector(j) = rhoLevel1Vectortemp1(j) + rhoLevel1Vectortemp2(j);
    hbar.iterate(1);
end
close(hbar); EstimateTime = toc;

Trix =zeros(3,size(elementsFEM,1)); Triy =zeros(3,size(elementsFEM,1)); Tric =zeros(1,size(elementsFEM,1)); Trictemp1=Tric; Trictemp2=Tric; 
TriStrain =zeros(4,size(elementsFEM,1));
for tempj = 1:size(elementsFEM,1)
    Trix(1:3,tempj) = coordinatesFEM(elementsFEM(tempj,1:3),1);
    Triy(1:3,tempj) = coordinatesFEM(elementsFEM(tempj,1:3),2);
    Tric(tempj) = rhoLevel1Vector(tempj); Trictemp1(tempj) = rhoLevel1Vectortemp1(tempj); Trictemp2(tempj) = rhoLevel1Vectortemp2(tempj);
    TriStrain(1,tempj) = CSTStrain(tempj,1); TriStrain(2,tempj) = CSTStrain(tempj,2);
    TriStrain(3,tempj) = CSTStrain(tempj,3); TriStrain(4,tempj) = CSTStrain(tempj,4);
end
close all;
figure; patch(Trix,Triy,Tric,'edgecolor','none' ); c=colorbar; title('Ele total err estimator','fontweight','normal'); set(gca,'fontsize',16); axis tight; axis equal;  set(gca,'XTick',[]);
figure; patch(Trix,Triy,Trictemp1,'edgecolor','none');c=colorbar; title('Ele interior err estimator','fontweight','normal'); set(gca,'fontsize',16); axis tight; axis equal;  set(gca,'XTick',[]);
figure; patch(Trix,Triy,Trictemp2,'edgecolor','none');c=colorbar; title('Ele jump err estimator','fontweight','normal'); set(gca,'fontsize',16); axis tight; axis equal;  set(gca,'XTick',[]);


figure(1); savefig(['fig_L',num2str(LevelNo),'_err_all.fig']);
figure(2); savefig(['fig_L',num2str(LevelNo),'_err_i.fig']);
figure(3); savefig(['fig_L',num2str(LevelNo),'_err_j.fig']);


% ======= Compute SSD Error and Penalty Error =======
[ErrSSD] = funSSDErrTri(coordinatesFEM,elementsFEM,fNormalized,gNormalized,ULevel1,0,0);
sum(ErrSSD(:))
% % ------- Plot SSD ---------
% figure; surf(reshape(ErrSSD,M-1,N-1)'); view(2); axis tight; axis equal; colormap jet; box on;
% title('SSD correlation function','fontweight','normal','Interpreter','latex'); set(gca,'fontsize',18); 
% xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex'); set(gcf,'color','w');
% a = gca; a.TickLabelInterpreter = 'latex'; b = colorbar; b.TickLabelInterpreter = 'latex';

figure; patch(Trix,Triy,ErrSSD,'edgecolor','none');c=colorbar; title('Ele err estimator AL1','fontweight','normal'); set(gca,'fontsize',16); axis tight; axis equal;  set(gca,'XTick',[]);
savefig(['fig_L',num2str(LevelNo),'_err_AL1.fig']); % JY!!!


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DorflerSatisfied = 1; 
while DorflerSatisfied == 1
    fprintf('--- Dorfler strategy to mark elements to refine? --- \n')
    fprintf('Input a real number between 0 and 1 here: ');
    prompt = 'Input here: '; thetaDorfler = input(prompt);
    
    
    
    %% ================ Mark: Level1 ==================
    % ======= Dorfler's Strategy to Mark elements to refine =======
    tic; [rhoLevel1VectorSorted,rhoLevel1VectorSortIndex] = sort(rhoLevel1Vector,'descend');
    thetaDorflerList(LevelNo) = thetaDorfler; rhoLevel1SumM = 0; rhoLevel1SumAll = sum(rhoLevel1Vector); rhoLevel1SumMIndex = 0;
    while rhoLevel1SumM < 0.999999*thetaDorfler*rhoLevel1SumAll
        rhoLevel1SumMIndex = rhoLevel1SumMIndex + 1;
        rhoLevel1SumM = rhoLevel1SumM + rhoLevel1VectorSorted(rhoLevel1SumMIndex);
    end
    refinementPatch = rhoLevel1VectorSortIndex(1:rhoLevel1SumMIndex);
    if CrackOrNot > 0
        refinementPatch = union(refinementPatch,[EleCrackTop(end)+1:1:EleCrackBottom(1)-1]');
    end
    Trix = zeros(3,size(elementsFEM,1)); Triy = zeros(3,size(elementsFEM,1)); Tric = zeros(1,size(elementsFEM,1));
    for j = 1:size(elementsFEM,1)
        Trix(1:3,j) = coordinatesFEM(elementsFEM(j,1:3),1);
        Triy(1:3,j) = coordinatesFEM(elementsFEM(j,1:3),2);
    end
    for j = 1:size(refinementPatch,1),Tric(refinementPatch(j)) = 1;end
    MarkTime = toc;
    figure; patch(Trix,Triy,Tric,'EdgeColor','none'); title('Marked ele (yellow)','fontweight','normal'); set(gca,'fontsize',16); axis tight; axis equal;
 
    savefig(['fig_L',num2str(LevelNo+1),'_mark.fig']); % JY!!!
    
    
    %% ============== Refine: Level1 ==============
    if thetaDorfler>0 && size(coordinatesFEM,1)>1   
        
        tic; [coordinatesFEMLevel2,elementsFEMLevel2,dirichletLevel2,neumannLevel2] = ...
            TrefineNVBwBCNormal(coordinatesFEM,elementsFEM,dirichlet,neumann,refinementPatch);
        RefineTime = toc;
     
        % Assign element generation #
        eleGenerationLevel2 = ones(size(elementsFEMLevel2,1),1);
        for tempj=1:size(elementsFEMLevel2,1)
            vol0 = 0.25*winstepsize^2;
            voltemp = 0.5*det([ coordinatesFEMLevel2(elementsFEMLevel2(tempj,1),:),1;
                                coordinatesFEMLevel2(elementsFEMLevel2(tempj,2),:),1;
                                coordinatesFEMLevel2(elementsFEMLevel2(tempj,3),:),1  ]);
            eleGenerationLevel2(tempj) = 1-log2( voltemp/vol0 );
        end
        
        % Initialize variables F0_Level2 and U0_Level2
        F_dispu = scatteredInterpolant( coordinatesFEM(:,1),coordinatesFEM(:,2),ULevel1(1:2:end) );
        F_dispv = scatteredInterpolant( coordinatesFEM(:,1),coordinatesFEM(:,2),ULevel1(2:2:end) );
        F_F11 = scatteredInterpolant( coordinatesFEM(:,1),coordinatesFEM(:,2),FLevel1(1:4:end) );
        F_F21 = scatteredInterpolant( coordinatesFEM(:,1),coordinatesFEM(:,2),FLevel1(2:4:end) );
        F_F12 = scatteredInterpolant( coordinatesFEM(:,1),coordinatesFEM(:,2),FLevel1(3:4:end) );
        F_F22 = scatteredInterpolant( coordinatesFEM(:,1),coordinatesFEM(:,2),FLevel1(4:4:end) );
        
        U0Level2 = 0*coordinatesFEMLevel2(:); F0Level2=0*[U0Level2(:);U0Level2(:)];
        temp = F_dispu(coordinatesFEMLevel2(:,1),coordinatesFEMLevel2(:,2)); U0Level2(1:2:end)=temp(:);
        temp = F_dispv(coordinatesFEMLevel2(:,1),coordinatesFEMLevel2(:,2)); U0Level2(2:2:end)=temp(:);
        temp = F_F11(coordinatesFEMLevel2(:,1),coordinatesFEMLevel2(:,2)); F0Level2(1:4:end)=temp(:);
        temp = F_F21(coordinatesFEMLevel2(:,1),coordinatesFEMLevel2(:,2)); F0Level2(2:4:end)=temp(:);
        temp = F_F12(coordinatesFEMLevel2(:,1),coordinatesFEMLevel2(:,2)); F0Level2(3:4:end)=temp(:);
        temp = F_F22(coordinatesFEMLevel2(:,1),coordinatesFEMLevel2(:,2)); F0Level2(4:4:end)=temp(:);
         
    else
        
        tic; refinedEleIDListLevel1 = 0;
        coordinatesFEMLevel2 = coordinatesFEM; elementsFEMLevel2 = elementsFEM; U0Level2 = ULevel1; F0Level2 = FLevel1;
        dirichletLevel2 = unique(dirichlet(:)); neumannLevel2 = neumann; eleGenerationLevel2 = eleGeneration;
        h=waitbar(0,'Wait for mesh refinement.');
        %tempNeuMat = [];
        for j =  1:size(refinementPatch,1)
            checkAlreadyRefinedOrNot = ismember(refinementPatch(j),refinedEleIDListLevel1);
            if checkAlreadyRefinedOrNot == 0
                FinishRefineRecursive = 1;
                [coordinatesFEMLevel2,elementsFEMLevel2,U0Level2,F0Level2,refinedEleIDListNew,dirichletNew,neumannNew] = ...
                    refineRecursiveTri(coordinatesFEMLevel2,elementsFEMLevel2(:,1:3),U0Level2,F0Level2,refinementPatch(j), refinementPatch, FinishRefineRecursive);

                refinedEleIDListLevel1 = unique([refinedEleIDListLevel1;refinedEleIDListNew]);
                dirichletLevel2 = unique([dirichletLevel2;dirichletNew']);
                % this step is wrong. newmann size is L*4: 
                if size(neumannNew,1)>1 
                    %tempNeuMat = [tempNeuMat; refinementPatch(j)];
                    [rowtemp1,~] = find( (neumannLevel2(:,1) == neumannNew(2,1))  );
                    [rowtemp2,~] = find( (neumannLevel2(:,2) == neumannNew(2,2))  );
                    [rowtemp3,~] = find( (neumannLevel2(:,2) == neumannNew(2,1))  );
                    [rowtemp4,~] = find( (neumannLevel2(:,1) == neumannNew(2,2))  );
                    rowtemp = union(intersect(rowtemp1,rowtemp2), intersect(rowtemp3,rowtemp4))  ; 
                    if neumannLevel2(rowtemp,1) ~= neumannNew(2,3),
                        LengthOfNeumanntemp = size(neumannLevel2,1);
                        neumannLevel2(LengthOfNeumanntemp+1,1:4) = [neumannLevel2(rowtemp,1),neumannNew(2,3),neumannLevel2(rowtemp,3:4)];
                        neumannLevel2(rowtemp,1) = [neumannNew(2,3)];
                    end
                end
                % neumannLevel2 = unique([neumannLevel2;neumannNew']);
            end
            waitbar(j/size(refinementPatch,1));
        end
        RefineTime = toc; close(h); dirichletLevel2 = dirichletLevel2(2:end);

        if CrackOrNot > 0
            [CoordCrackTopNextIter,~,CoordCrackBottomNextIter,~] = funCrackTopOrBottom(coordinatesFEMLevel2,elementsFEMLevel2,CrackPath1,CrackPath2,CrackTip);
            F0Level2 = funSmoothStrainCrack(F0Level2,coordinatesFEMLevel2,elementsFEMLevel2,winstepsize,CoordCrackTopNextIter,CoordCrackBottomNextIter,LevelNo+1,0,0);
        end

        dirichletLevel2 = [];
        
    end % END of mesh refinement

    % ====== Show refined mesh ======
    Trix = zeros(3,size(elementsFEMLevel2,1)); Triy = zeros(3,size(elementsFEMLevel2,1)); Tric = zeros(1,size(elementsFEMLevel2,1));
    for j = 1:size(elementsFEMLevel2,1)
        Trix(1:3,j) = coordinatesFEMLevel2(elementsFEMLevel2(j,1:3),1);
        Triy(1:3,j) = coordinatesFEMLevel2(elementsFEMLevel2(j,1:3),2);
        Tric(j) =  1;
    end
    figure; patch(Trix,Triy,Tric);  axis equal; axis off; title('Level2 Mesh','fontweight','normal'); set(gca,'fontsize',16);

    
    savefig(['fig_L',num2str(LevelNo+1),'_refine.fig']); % JY!!!
    
    
    fprintf('--- Are you satisfied with current Dorfler refinement number (0-yes; 1-no)? ---  \n')
    prompt = 'Input here: ';
    DorflerSatisfied = input(prompt);
    if DorflerSatisfied == 0
        fprintf('--- Mesh refinement done! ---  \n')
    end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== using fft to find integer research of new nodes ======
% PhiLevel2 = zeros((size(coordinatesFEMLevel2,1)-size(coordinatesFEM,1)),1);
% for tempi = 1:(size(coordinatesFEMLevel2,1)-size(coordinatesFEM,1))
%     tempindex = size(coordinatesFEM,1) + tempi ;
%     coordxtemp = coordinatesFEMLevel2(tempindex,1); coordytemp = coordinatesFEMLevel2(tempindex,2);
%     C = fNormalized(coordxtemp-winsize/2:coordxtemp+winsize/2, coordytemp-winsize/2:coordytemp+winsize/2);
%     D = gNormalized(coordxtemp-winsize/2-tempSizeOfSearchRegion:coordxtemp+winsize/2+tempSizeOfSearchRegion, ...
%         coordytemp-winsize/2-tempSizeOfSearchRegion:coordytemp+winsize/2+tempSizeOfSearchRegion);
%     XCORRF2OFCD0 = normxcorr2(C,D);
%     [v1temp,u1temp,max_f] = findpeak(XCORRF2OFCD0(winsize:end-winsize+1, winsize:end-winsize+1),1);
%     zero_disp = ceil(size(XCORRF2OFCD0(winsize:end-winsize+1,winsize:end-winsize+1))/2);
%     
%     PhiLevel2(tempi)        = max_f;
%     if max_f > 0.9
%     U0Level2(2*tempindex-1) = u1temp-zero_disp(1);
%     U0Level2(2*tempindex)   = v1temp-zero_disp(2);
%     end 
% end

% ====== inpaint nans for new nodes using gridfit ======
% Coordxnodes = unique(coordinatesFEMLevel2(:,1)); Coordynodes = unique(coordinatesFEMLevel2(:,2));
% [u1temp] = gridfit(coordinatesFEMLevel2(:,1), coordinatesFEMLevel2(:,2), U0Level2(1:2:end), Coordxnodes,Coordynodes); u1temp = u1temp';
% [v1temp] = gridfit(coordinatesFEMLevel2(:,1), coordinatesFEMLevel2(:,2), U0Level2(2:2:end), Coordxnodes,Coordynodes); v1temp = v1temp';
% for tempi = 1:(size(coordinatesFEMLevel2,1)-size(coordinatesFEM,1))
%     tempindex = size(coordinatesFEM,1) + tempi ;
%     [row1,col1] = find(Coordxnodes==coordinatesFEMLevel2(tempindex,1));
%     [row2,col2] = find(Coordynodes==coordinatesFEMLevel2(tempindex,2));
%     U0Level2(2*tempindex-1) = u1temp(row1,row2);
%     U0Level2(2*tempindex)   = v1temp(row1,row2);
% end

% ====== Plot refined mesh ======
%figure; show(elementsFEMLevel2,[],coordinatesFEMLevel2,U0Level2(1:2:end));
%title('x-Disp U','fontweight','normal'); set(gca,'fontsize',16); view(2); axis tight; axis equal; colorbar; set(gca,'XTick',[]);
% figure; show(elementsFEMLevel2,[],coordinatesFEMLevel2,U0Level2(2:2:end));
% title('y-Disp V','fontweight','normal'); set(gca,'fontsize',16); view(2); axis tight; axis equal; colorbar; set(gca,'XTick',[]);
% figure; show(elementsFEMLevel2,[],coordinatesFEMLevel2,F0Level2(1:4:end));
% title('F11','fontweight','normal'); set(gca,'fontsize',16); view(2); axis tight; axis equal; colorbar; set(gca,'XTick',[]);
% figure; show(elementsFEMLevel2,[],coordinatesFEMLevel2,F0Level2(4:4:end));
% title('F22','fontweight','normal'); set(gca,'fontsize',16); view(2); axis tight; axis equal; colorbar; set(gca,'XTick',[]);
% 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update winsize w/ {U0Level2, F0Level2}
LevelNo = LevelNo + 1; %winsizeList = winsize*ones(size(coordinatesFEMLevel2,1),2);
% 
% Ftemp1 = F0Level2(1:2:end); Ftemp2 = F0Level2(2:2:end);
% % DFtemp1 = {F11,x F12,x F11,y F12,y} = {u,xx u,yx u,xy u,yy};
% [DFtemp1,~,~] = Global_CSTStrain(coordinatesFEMLevel2,elementsFEMLevel2,Ftemp1);
% 
% %Plotstrain_showTri(DFtemp1,coordinatesFEMLevel2,elementsFEMLevel2(:,1:3));
% 
% winsizex_ub1 = abs(2*F0Level2(1:4:end)./DFtemp1(1:4:end));
% winsizex_ub2 = abs(2*F0Level2(3:4:end)./DFtemp1(3:4:end));
% winsizey_ub1 = abs(2*F0Level2(1:4:end)./DFtemp1(2:4:end));
% winsizey_ub2 = abs(2*F0Level2(3:4:end)./DFtemp1(4:4:end));
% 
% % DFtemp2 = {F21,x F22,x F21,y F22,y} = {v,xx v,yx v,xy v,yy};
% [DFtemp2,~,~] = Global_CSTStrain(coordinatesFEMLevel2,elementsFEMLevel2,Ftemp2);
% 
% winsizex_ub3 = abs(2*F0Level2(2:4:end)./DFtemp2(1:4:end));
% winsizex_ub4 = abs(2*F0Level2(4:4:end)./DFtemp2(3:4:end));
% winsizey_ub3 = abs(2*F0Level2(2:4:end)./DFtemp2(2:4:end));
% winsizey_ub4 = abs(2*F0Level2(4:4:end)./DFtemp2(4:4:end));
% 
% winsizex_ub = round(min([winsizex_ub1,winsizex_ub2,winsizex_ub3,winsizex_ub4,winsize*ones(length(winsizex_ub1),1)],[],2));
% winsizexList = max([winsizex_ub,16*ones(length(winsizex_ub1),1)],[],2);
% winsizey_ub = round(min([winsizey_ub1,winsizey_ub2,winsizey_ub3,winsizey_ub4,winsize*ones(length(winsizey_ub1),1)],[],2));
% winsizeyList = max([winsizey_ub,16*ones(length(winsizey_ub1),1)],[],2);
% winsizeList = 2*ceil([winsizexList,winsizeyList]/2);
% winsizeListT = winsizeList';
% 
% Plotdisp_showTri(winsizeListT(:),coordinatesFEMLevel2,elementsFEMLevel2);
% fign = get(gcf,'Number'); figure(fign), title('Subset $y$-direction size (pixels)','FontWeight','Normal','Interpreter','latex')
% figure(fign-1), title('Subset $x$-direction size (pixels)','FontWeight','Normal','Interpreter','latex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% ============== Level2 ==============        
PassCrackOrNotIter = zeros(size(elementsFEMLevel2,1),1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ULevel2, normOfWLevel2, TimeICGNLevel2] = funGlobalICGNTri(coordinatesFEMLevel2,elementsFEMLevel2,Df,fNormalized,gNormalized,U0Level2,alphavar,tol,ClusterNo);
[FLevel2] = funGlobal_CSTStrain(coordinatesFEMLevel2,elementsFEMLevel2,ULevel2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[USubpb1,FSubpb1] = funRemoveOutliersTri_adapt(coordinatesFEMLevel2,elementsFEMLevel2,USubpb1,FSubpb1,0);
% close all; Plotdisp_showTri(USubpb1,coordinatesFEMLevel2,elementsFEMLevel2);
% Plotstrain_showTri(FSubpb1,coordinatesFEMLevel2,elementsFEMLevel2);
%save(['Subpb1_step',num2str(ALSolveStep)],'USubpb1','FSubpb1');

if LevelNo < 6
    if CrackOrNot > 0
%         % ------- Find elements located on the top/bottom of crack path -------                        
%         [CoordCrackTop,EleCrackTop,CoordCrackBottom,EleCrackBottom] = funCrackTopOrBottom(coordinatesFEMLevel2,elementsFEMLevel2,CrackPath1,CrackPath2,CrackTip);
%         % ------- Before computing strain, we smooth the displacement field -------
%         USubpb2 = funSmoothDispCrack(ULocal,coordinatesFEMLevel2,elementsFEMLevel2,winstepsize,CoordCrackTop,CoordCrackBottom,LevelNo,DispFilterSize,DispFilterStd);
%         % -------- Plot disp -------- 
%         Plotdisp_show(USubpb2,coordinatesFEMLevel2,elementsFEMLevel2([EleCrackTop;EleCrackBottom],1:4)); 
%         % -------- Compute strain field -------- 
%         FLocal = ComputeCrackStrain(USubpb2,FSubpb2,coordinatesFEMLevel2,elementsFEMLevel2,EleCrackTop,EleCrackBottom,EnrHAndTipEleIndex,EnrTipEleIndex,0,CrackPath1,CrackPath2,CrackTip,0);
%         FSubpb2 = funSmoothStrainCrack(FLocal,coordinatesFEMLevel2,elementsFEMLevel2,winstepsize,CoordCrackTop,CoordCrackBottom,LevelNo,StrainFilterSize,StrainFilterStd);
%         Plotstrain_show(FSubpb2,coordinatesFEMLevel2,elementsFEMLevel2([EleCrackTop;EleCrackBottom],1:4));
    else
%         USubpb1 = funSmoothDisp_adapt(USubpb1,coordinatesFEMLevel2,elementsFEMLevel2,winstepsize,LevelNo,DispFilterSize,DispFilterStd);
%         Plotdisp_show(USubpb1,coordinatesFEMLevel2,elementsFEMLevel2(:,1:4));
%         FSubpb1 = funSmoothStrain_adapt(FSubpb1,coordinatesFEMLevel2,elementsFEMLevel2,winstepsize,LevelNo,StrainFilterSize,StrainFilterStd);
%         Plotstrain_show(FSubpb1,coordinatesFEMLevel2,elementsFEMLevel2(:,1:4));
    end
end

%[FSubpb12,~,~] = Global_NodalStrainAvg(coordinatesFEMLevel2,elementsFEMLevel2,USubpb1,2,...
%                    winstepsize,LevelNo,EnrHAndTipEleIndex,EnrTipEleIndex,CrackOrNot,CrackPath1,CrackPath2,CrackTip,CrackTipOrNot);
%        Plotstrain_show(FSubpb12,coordinatesFEMLevel2,elementsFEMLevel2);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if CrackOrNot > 0
%     [USubpb2,~] = ComputeCrackDisp(USubpb2temp,x0,y0,coordinatesFEMLevel2,elementsFEMLevel2,EnrTipEleIndex,CrackOrNot,CrackPath1,CrackPath2,CrackTip,CrackTipOrNot);
%     USubpb2 = full(USubpb2); % close all; Plotdisp_show(USubpb2,coordinatesFEMLevel2,elementsFEMLevel2);
% else, USubpb2 = full(USubpb2temp); % close all; Plotdisp_show(USubpb2,coordinatesFEMLevel2,elementsFEMLevel2);
% end
% ======== Smooth the displacement field =======
try
  temp = DispFilterSize; temp = DispFilterStd; temp = StrainFilterSize; temp = StrainFilterStd;
catch
 DispFilterSize = 0; DispFilterStd = 0; StrainFilterSize = 0; StrainFilterStd = 0;
end
% USubpb2 = funSmoothDisp_adapt(USubpb2,M,N,x,y);
%[FSubpb2,~,~] = Global_CSTStrain(coordinatesFEMLevel2,elementsFEMLevel2,USubpb2 );

% [FSubpb2,StrainGaussPt,CoordGaussPt] = Global_NodalStrainAvg(coordinatesFEMLevel2,elementsFEMLevel2,USubpb1,2,...
%     winstepsize,LevelNo,EnrHAndTipEleIndex,EnrTipEleIndex,CrackOrNot,CrackPath1,CrackPath2,CrackTip,CrackTipOrNot);
% tritemp = delaunay(CoordGaussPt(:,1),CoordGaussPt(:,2));
% figure,trisurf(tritemp,CoordGaussPt(:,1),CoordGaussPt(:,2),StrainGaussPt(:,1));

% ------- Smooth strain field --------
%FSubpb2 = funSmoothStrain_adapt(FSubpb2,coordinatesFEMLevel2,elementsFEMLevel2,winstepsize/2,LevelNo,StrainFilterSize,StrainFilterStd);
  

% Plotdisp_show(USubpb1,coordinatesFEMLevel2,elementsFEMLevel2(:,1:4));
% Plotstrain_show(FSubpb1,coordinatesFEMLevel2,elementsFEMLevel2(:,1:4));
% Plotdisp_show(USubpb2,coordinatesFEMLevel2,elementsFEMLevel2(:,1:4));
% Plotstrain_show(FSubpb2,coordinatesFEMLevel2,elementsFEMLevel2(:,1:4));
  
if LevelNo < 6
    if CrackOrNot > 0
%         % ------- Find elements located on the top/bottom of crack path -------                        
%         [CoordCrackTop,EleCrackTop,CoordCrackBottom,EleCrackBottom] = funCrackTopOrBottom(coordinatesFEMLevel2,elementsFEMLevel2,CrackPath1,CrackPath2,CrackTip);
%         % ------- Before computing strain, we smooth the displacement field -------
%         USubpb2 = funSmoothDispCrack(USubpb2,coordinatesFEMLevel2,elementsFEMLevel2,winstepsize,CoordCrackTop,CoordCrackBottom,LevelNo,DispFilterSize,DispFilterStd);
%         % -------- Plot disp -------- 
%         Plotdisp_show(USubpb2,coordinatesFEMLevel2,elementsFEMLevel2([EleCrackTop;EleCrackBottom],1:4)); 
%         % -------- Compute strain field -------- 
%         FLocal = ComputeCrackStrain(USubpb2,FSubpb2,coordinatesFEMLevel2,elementsFEMLevel2,EleCrackTop,EleCrackBottom,EnrHAndTipEleIndex,EnrTipEleIndex,0,CrackPath1,CrackPath2,CrackTip,0);
%         FSubpb2 = funSmoothStrainCrack(FLocal,coordinatesFEMLevel2,elementsFEMLevel2,winstepsize,CoordCrackTop,CoordCrackBottom,LevelNo,StrainFilterSize,StrainFilterStd);
%         Plotstrain_show(FSubpb2,coordinatesFEMLevel2,elementsFEMLevel2([EleCrackTop;EleCrackBottom],1:4));
    else
        close all;
        % USubpb2 = funSmoothDisp_adapt(USubpb2,coordinatesFEMLevel2,elementsFEMLevel2,winstepsize,LevelNo,DispFilterSize,DispFilterStd);
         Plotdisp_showTri(ULevel2,coordinatesFEMLevel2,elementsFEMLevel2(:,1:3));
        % FSubpb2 = funSmoothStrain_adapt(FSubpb2,coordinatesFEMLevel2,elementsFEMLevel2,winstepsize,LevelNo,StrainFilterSize,StrainFilterStd);
          Plotstrain_showTri(FLevel2,coordinatesFEMLevel2,elementsFEMLevel2(:,1:3));
          
          figure(1); savefig(['fig_L',num2str(LevelNo),'_dispu.fig']);
          figure(2); savefig(['fig_L',num2str(LevelNo),'_dispv.fig']);
          figure(3); savefig(['fig_L',num2str(LevelNo),'_e11.fig']);
          figure(4); savefig(['fig_L',num2str(LevelNo),'_e22.fig']);
          figure(5); savefig(['fig_L',num2str(LevelNo),'_e12.fig']);
          figure(6); savefig(['fig_L',num2str(LevelNo),'_eshear.fig']);
          
    end
end
disp(['--- Done! step ',num2str(LevelNo),' subproblem 2 ---']);
 