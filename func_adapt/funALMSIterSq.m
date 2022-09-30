%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that do adaptive ALDIC 
% in the oder of: (Input: (LevelNo-1)Solve) -> Estimate -> Mark -> Refine -> Solve(Output)
% 
% In this code, variables at (LevelNo-1) are all called "Level1", while
% variables at (LevelNo) are called "Level2".
%
% Author: Jin Yang, Email: jyang526@wisc.edu
% Date: 2019 June.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [UIter,FIter,udualIter,vdualIter,USubpb2tempIter,ConvItPerEle,betamuIter,...
    coordinatesFEMLevel2,elementsFEMLevel2,eleGenerationLevel2,irregularEdgeLevel2,dirichletLevel2,neumannLevel2,...
    thetaDorfler,ALSub1Time,ALSub2Time,EstimateTime,MarkTime,RefineTime,...
    rhoLevel1Vector,rhoLevel1Vectortemp1,rhoLevel1Vectortemp2,EnrHAndTipEleIndex,EnrTipEleIndex] ...
    = funALMSIterSq(ULevel1,FLevel1,udualLevel1,vdualLevel1,USubpb2tempLevel1,...
    coordinatesFEM,elementsFEM,winstepsize,ClusterNo,...
    eleGeneration,irregularEdge,dirichlet,neumann,LevelNo,fNormalized,gNormalized,winsize,...
    CrackOrNot,CrackPath1,CrackPath2,CrackTip,CrackTipOrNot,EnrHAndTipEleIndex,EnrTipEleIndex,...
    DispFilterSize,DispFilterStd,StrainFilterSize,StrainFilterStd,tol,tol2,betamuLevel1)

% ULevel1 = UIter; FLevel1 = FIter; udualLevel1 = udualIter; vdualLevel1 = vdualIter; USubpb2tempLevel1 = USubpb2tempIter;
% elementsFEM = elementsFEMIter;  coordinatesFEM = coordinatesFEMIter;
% eleGeneration = eleGenerationIter; dirichlet = dirichletIter; neumann = neumannIter; irregularEdge=zeros(1,3);
% fNormalized = ImgNormalized{1}; gNormalized = ImgNormalized{2};
% betamuLevel1 = betamuIter;
warning('off');
betavar = betamuLevel1.betavar; muvar = betamuLevel1.mu; alphavar = betamuLevel1.alphavar;
imgPyramidUnit = 1;

%% ============== Estimate: Level1 ==============
% ------ Compute rho to find where to refine ------
rhoLevel1Vector = zeros(size(elementsFEM,1),1); rhoLevel1Vectortemp1 = rhoLevel1Vector; rhoLevel1Vectortemp2 = rhoLevel1Vector;
ULevel1temp = [ULevel1;0;0]; hbar = parfor_progressbar(size(elementsFEM,1),'Please wait for a posterior error estimate!'); tic

parfor j = 1:size(elementsFEM,1)
    % eleNeighborIndexAndEdge = findEleNeighbors(elementsFEMIter, coordinatesFEMIter, 23 ,eleGenerationIter, winstepsize);
    point1x = coordinatesFEM(elementsFEM(j,1),1); point3x = coordinatesFEM(elementsFEM(j,3),1);lengthOfElement = (point3x-point1x)/winstepsize;
    % ---- Compute error estimator ----
    rhoLevel1Vectortemp1(j) = (lengthOfElement)^2 * aPostErrEstInt_ALSq(coordinatesFEM,elementsFEM,j,betavar,muvar,ULevel1,FLevel1,udualLevel1,vdualLevel1,imgPyramidUnit,CrackOrNot);
    rhoLevel1Vectortemp2(j) = (lengthOfElement) * (betavar+alphavar)^2 * aPostErrEstJump_ALSq(coordinatesFEM,elementsFEM,j,eleGeneration,winstepsize,USubpb2tempLevel1, ...
                               EnrHAndTipEleIndex,EnrTipEleIndex,CrackOrNot,CrackPath1,CrackPath2,CrackTip,CrackTipOrNot);
    rhoLevel1Vector(j) = rhoLevel1Vectortemp1(j) + rhoLevel1Vectortemp2(j);
    hbar.iterate(1);
end
close(hbar); EstimateTime = toc;
Sqx = zeros(4,size(elementsFEM,1)); Sqy = zeros(4,size(elementsFEM,1)); Sqc = zeros(1,size(elementsFEM,1)); Sqctemp1 = Sqc; Sqctemp2 = Sqc;
for tempj = 1:size(elementsFEM,1)
    Sqx(1:4,tempj) = coordinatesFEM(elementsFEM(tempj,1:4),1);Sqy(1:4,tempj) = coordinatesFEM(elementsFEM(tempj,1:4),2); 
    Sqc(tempj) = rhoLevel1Vector(tempj);Sqctemp1(tempj) = rhoLevel1Vectortemp1(tempj);Sqctemp2(tempj) = rhoLevel1Vectortemp2(tempj);
end
close all;
figure; patch(Sqx,Sqy,Sqc,'edgecolor','none' ); c=colorbar; title('Ele total err estimator','fontweight','normal'); set(gca,'fontsize',16); axis tight; axis equal;  set(gca,'XTick',[]);
figure; patch(Sqx,Sqy,Sqctemp1,'edgecolor','none');c=colorbar; title('Ele interior err estimator','fontweight','normal'); set(gca,'fontsize',16); axis tight; axis equal;  set(gca,'XTick',[]);
figure; patch(Sqx,Sqy,Sqctemp2,'edgecolor','none');c=colorbar; title('Ele jump err estimator','fontweight','normal'); set(gca,'fontsize',16); axis tight; axis equal;  set(gca,'XTick',[]);
        
figure(1); savefig(['fig_L',num2str(LevelNo),'_err_all.fig']);
figure(2); savefig(['fig_L',num2str(LevelNo),'_err_i.fig']);
figure(3); savefig(['fig_L',num2str(LevelNo),'_err_j.fig']);

% ======= Compute SSD Error and Penalty Error ======= % JY!!!
[ErrSSD] = funSSDErr(coordinatesFEM,elementsFEM,fNormalized,gNormalized,ULevel1,0,0);
sum(ErrSSD(:))
% ------- Plot SSD ---------
% figure; surf(reshape(ErrSSD,M-1,N-1)'); view(2); axis tight; axis equal; colormap jet; box on;
% title('SSD correlation function','fontweight','normal','Interpreter','latex'); set(gca,'fontsize',18); 
% xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex'); set(gcf,'color','w');
% a = gca; a.TickLabelInterpreter = 'latex'; b = colorbar; b.TickLabelInterpreter = 'latex';
 
figure; patch(Sqx,Sqy,ErrSSD,'edgecolor','none');c=colorbar; title('Ele err estimator AL1','fontweight','normal'); set(gca,'fontsize',16); axis tight; axis equal;  set(gca,'XTick',[]);
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
    rhoLevel1SumM = 0; rhoLevel1SumAll = sum(rhoLevel1Vector(:)); rhoLevel1SumMIndex = 0;
    while rhoLevel1SumM < 0.999999*thetaDorfler*rhoLevel1SumAll
        rhoLevel1SumMIndex = rhoLevel1SumMIndex + 1;
        rhoLevel1SumM = rhoLevel1SumM + rhoLevel1VectorSorted(rhoLevel1SumMIndex);
    end
    refinementPatch = rhoLevel1VectorSortIndex(1:rhoLevel1SumMIndex);
    if CrackOrNot > 0
        refinementPatch = union(refinementPatch,[EleCrackTop(end)+1:1:EleCrackBottom(1)-1]');
    end
    Sqx = zeros(4,size(elementsFEM,1)); Sqy = zeros(4,size(elementsFEM,1)); Sqc = zeros(1,size(elementsFEM,1));
    for j = 1:size(elementsFEM,1); Sqx(1:4,j) = coordinatesFEM(elementsFEM(j,1:4),1); Sqy(1:4,j) = coordinatesFEM(elementsFEM(j,1:4),2); end
    for j = 1:size(refinementPatch,1); Sqc(refinementPatch(j)) = 1; end
    MarkTime = toc;
    figure; patch(Sqx,Sqy,Sqc,'EdgeColor','none'); caxis([0,1]); colorbar;
    title('Marked ele (yellow)','fontweight','normal'); set(gca,'fontsize',16); axis tight; axis equal;  

    savefig(['fig_L',num2str(LevelNo+1),'_mark.fig']); % JY!!!
 
    
    %% ============== Refine: Level1 ==============
    if thetaDorfler>0 && size(coordinatesFEM,1)>1e3  
    
        tic; [coordinatesFEMLevel2,elementsFEMLevel2,irregularEdgeLevel2,dirichletLevel2,neumannLevel2] = ...
            QrefineRwBCNormal(coordinatesFEM,elementsFEM(:,1:4),irregularEdge,dirichlet,neumann(:,1:4),refinementPatch);
        elementsFEMLevel2 = [elementsFEMLevel2,zeros(size(elementsFEMLevel2,1),4)];
        
        % Assign element generation #
        eleGenerationLevel2 = ones(size(elementsFEMLevel2,1),1);
        for tempj=1:size(elementsFEMLevel2,1)
            eleGenerationLevel2(tempj) = 1-log2(norm(coordinatesFEMLevel2(elementsFEMLevel2(tempj,1),:) - coordinatesFEMLevel2(elementsFEMLevel2(tempj,2),:))/winstepsize);
        end
         
        % Reorder element nodes
        for tempj=1:length(eleGenerationLevel2)
            [~,row2] = max(  sum(coordinatesFEMLevel2(elementsFEMLevel2(tempj,1:4)',:)')  );
            if row2==1
                elementsFEMLevel2(tempj,1:4) = elementsFEMLevel2(tempj,[3,4,1,2]);
            elseif row2==2
                elementsFEMLevel2(tempj,1:4) = elementsFEMLevel2(tempj,[4,1,2,3]);
            elseif row2==4
                elementsFEMLevel2(tempj,1:4) = elementsFEMLevel2(tempj,[2,3,4,1]);
            end
        end
        
        % For hanging nodes in FE-mesh
        for tempj=1:size(irregularEdgeLevel2,1)
            [Lia, Locb] = ismember( [irregularEdgeLevel2(tempj,1:2)], elementsFEMLevel2(:,1:2), 'rows' );
            if Lia>0, elementsFEMLevel2(Locb,8)=irregularEdgeLevel2(tempj,3); end
            [Lia, Locb] = ismember( [irregularEdgeLevel2(tempj,1:2)], elementsFEMLevel2(:,2:3), 'rows' );
            if Lia>0, elementsFEMLevel2(Locb,5)=irregularEdgeLevel2(tempj,3); end
            [Lia, Locb] = ismember( [irregularEdgeLevel2(tempj,1:2)], elementsFEMLevel2(:,3:4), 'rows' );
            if Lia>0, elementsFEMLevel2(Locb,6)=irregularEdgeLevel2(tempj,3); end
            [Lia, Locb] = ismember( [irregularEdgeLevel2(tempj,1:2)], elementsFEMLevel2(:,[4,1]), 'rows' );
            if Lia>0, elementsFEMLevel2(Locb,7)=irregularEdgeLevel2(tempj,3); end
            
            [Lia, Locb] = ismember( [irregularEdgeLevel2(tempj,1:2)], elementsFEMLevel2(:,[2,1]), 'rows' );
            if Lia>0, elementsFEMLevel2(Locb,8)=irregularEdgeLevel2(tempj,3); end
            [Lia, Locb] = ismember( [irregularEdgeLevel2(tempj,1:2)], elementsFEMLevel2(:,[3,2]), 'rows' );
            if Lia>0, elementsFEMLevel2(Locb,5)=irregularEdgeLevel2(tempj,3); end
            [Lia, Locb] = ismember( [irregularEdgeLevel2(tempj,1:2)], elementsFEMLevel2(:,[4,3]), 'rows' );
            if Lia>0, elementsFEMLevel2(Locb,6)=irregularEdgeLevel2(tempj,3); end
            [Lia, Locb] = ismember( [irregularEdgeLevel2(tempj,1:2)], elementsFEMLevel2(:,[1,4]), 'rows' );
            if Lia>0, elementsFEMLevel2(Locb,7)=irregularEdgeLevel2(tempj,3); end
        end
        RefineTime = toc;
        
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
                [coordinatesFEMLevel2,elementsFEMLevel2,U0Level2,F0Level2,refinedEleIDListNew,dirichletNew,neumannNew,eleGenerationLevel2] = ...
                    refineRecursiveSq(coordinatesFEMLevel2,elementsFEMLevel2,eleGenerationLevel2,U0Level2,F0Level2,refinementPatch(j),winstepsize,FinishRefineRecursive,...
                    CrackOrNot,CrackPath1,CrackPath2,CrackTip,CrackTipOrNot);
                
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
                    if neumannLevel2(rowtemp,1) ~= neumannNew(2,3)
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
        
        % Update irregularEdgeLevel2
        irregularEdgeLevel2 = [0,0,0];
        [row,~]=find(elementsFEMLevel2(:,5)>0);
        for tempi=1:length(row)
            irregularEdgeLevel2 = [irregularEdgeLevel2; elementsFEMLevel2(row(tempi),2:3),  elementsFEMLevel2(row(tempi),5) ];
        end
        [row,~]=find(elementsFEMLevel2(:,6)>0);
        for tempi=1:length(row)
            irregularEdgeLevel2 = [irregularEdgeLevel2; elementsFEMLevel2(row(tempi),3:4),  elementsFEMLevel2(row(tempi),6) ];
        end
        [row,~]=find(elementsFEMLevel2(:,7)>0);
        for tempi=1:length(row)
            irregularEdgeLevel2 = [irregularEdgeLevel2; elementsFEMLevel2(row(tempi),[4,1]),  elementsFEMLevel2(row(tempi),7) ];
        end
        [row,~]=find(elementsFEMLevel2(:,8)>0);
        for tempi=1:length(row)
            irregularEdgeLevel2 = [irregularEdgeLevel2; elementsFEMLevel2(row(tempi),1:2),  elementsFEMLevel2(row(tempi),8) ];
        end
        
        dirichletLevel2 = [];
         
    end % END of mesh refinement
    
    

    % ====== Show refined mesh ======
    Sqx = zeros(4,size(elementsFEMLevel2,1)); Sqy = zeros(4,size(elementsFEMLevel2,1)); Sqc = zeros(1,size(elementsFEMLevel2,1));
    for j = 1:size(elementsFEMLevel2,1)
        Sqx(1:4,j) = coordinatesFEMLevel2(elementsFEMLevel2(j,1:4),1);
        Sqy(1:4,j) = coordinatesFEMLevel2(elementsFEMLevel2(j,1:4),2); Sqc(j) =  1;
    end
    figure; patch(Sqx,Sqy,Sqc);  axis equal; axis off; title(['Level ',num2str(LevelNo+1),' Mesh'],'fontweight','normal'); set(gca,'fontsize',16);

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
% figure; show([], elementsFEMLevel2,coordinatesFEMLevel2,U0Level2(1:2:end));
% title('x-Disp U','fontweight','normal'); set(gca,'fontsize',16); view(2); axis tight; axis equal; colorbar; set(gca,'XTick',[]);
% figure; show([], elementsFEMLevel2,coordinatesFEMLevel2,U0Level2(2:2:end));
% title('y-Disp V','fontweight','normal'); set(gca,'fontsize',16); view(2); axis tight; axis equal; colorbar; set(gca,'XTick',[]);
% figure; show([], elementsFEMLevel2,coordinatesFEMLevel2,F0Level2(1:4:end));
% title('F11','fontweight','normal'); set(gca,'fontsize',16); view(2); axis tight; axis equal; colorbar; set(gca,'XTick',[]);
% figure; show([], elementsFEMLevel2,coordinatesFEMLevel2,F0Level2(4:4:end));
% title('F22','fontweight','normal'); set(gca,'fontsize',16); view(2); axis tight; axis equal; colorbar; set(gca,'XTick',[]);
% figure; show([], elementsFEMLevel2,coordinatesFEMLevel2,0.5*(F0Level2(2:4:end)+F0Level2(3:4:end)));
% title('0.5*(F12+F21)','fontweight','normal'); set(gca,'fontsize',16); view(2); axis tight; axis equal; colorbar; set(gca,'XTick',[]);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update winsize w/ {U0Level2, F0Level2}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LevelNo = LevelNo + 1;
% Ftemp1 = F0Level2(1:2:end); Ftemp2 = F0Level2(2:2:end);
% % DFtemp1 = {F11,x F12,x F11,y F12,y} = {u,xx u,yx u,xy u,yy};
% [DFtemp1,~,~] = Global_NodalStrainAvg(coordinatesFEMLevel2,elementsFEMLevel2,Ftemp1,2,...
%     winstepsize,LevelNo,EnrHAndTipEleIndex,EnrTipEleIndex,CrackOrNot,CrackPath1,CrackPath2,CrackTip,CrackTipOrNot);
% 
% Plotstrain_show(DFtemp1,coordinatesFEMLevel2,elementsFEMLevel2(:,1:4));
% 
% winsizex_ub1 = abs(2*F0Level2(1:4:end)./DFtemp1(1:4:end));
% winsizex_ub2 = abs(2*F0Level2(3:4:end)./DFtemp1(3:4:end));
% winsizey_ub1 = abs(2*F0Level2(1:4:end)./DFtemp1(2:4:end));
% winsizey_ub2 = abs(2*F0Level2(3:4:end)./DFtemp1(4:4:end));
% 
% % DFtemp2 = {F21,x F22,x F21,y F22,y} = {v,xx v,yx v,xy v,yy};
% [DFtemp2,~,~] = Global_NodalStrainAvg(coordinatesFEMLevel2,elementsFEMLevel2,Ftemp2,2,...
%     winstepsize,LevelNo,EnrHAndTipEleIndex,EnrTipEleIndex,CrackOrNot,CrackPath1,CrackPath2,CrackTip,CrackTipOrNot);
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
% Plotdisp_show(winsizeListT(:),coordinatesFEMLevel2,elementsFEMLevel2);
% fign = get(gcf,'Number'); figure(fign), title('Subset $y$-direction size (pixels)','FontWeight','Normal','Interpreter','latex')
% figure(fign-1), title('Subset $x$-direction size (pixels)','FontWeight','Normal','Interpreter','latex')
winsizeList=winsize;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% ============== Level2 ==============        
PassCrackOrNotIter = zeros(size(elementsFEMLevel2,1),1);
Df = funImgGradient(fNormalized,gNormalized); 

% ============== Subproblem 1: Level Iter ==============
% ------ Define penalty parameter ------
% mu = 1e-3; betavar = 0.1*winstepsize^2*mu; 
% mu = 0; betavar = 0; % tol = 1e-6; % When mu & betavar are 0s, ALDIC Sub1 is exactly the same as Local Subset DIC IC-GN 
ALSolveStep = 1; ALSub1Time = zeros(6,1); ALSub2Time = zeros(6,1); ConvItPerEle = zeros(size(coordinatesFEMLevel2,1),6);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------ Start Local DIC IC-GN iteration ------
[USubpb1,FSubpb1,HtempPar,ALSub1Timetemp,ConvItPerEletemp] = LocalICGN(U0Level2,coordinatesFEMLevel2,...
        Df,fNormalized,gNormalized,winsizeList,winstepsize/(2^(LevelNo-1)),tol,'GaussNewton',ClusterNo);
ALSub1Time(ALSolveStep) = ALSub1Timetemp; ConvItPerEle(:,ALSolveStep) = ConvItPerEletemp; clear ALSub1Timetemp ConvItPerEletemp; toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[USubpb1,FSubpb1] = funRemoveOutliers_adapt(coordinatesFEMLevel2,elementsFEMLevel2,USubpb1,FSubpb1,0);

% close all; Plotdisp_show(USubpb1,coordinatesFEMLevel2,elementsFEMLevel2);
% Plotstrain_show(FSubpb1,coordinatesFEMLevel2,elementsFEMLevel2);
save(['Subpb1_step',num2str(ALSolveStep)],'USubpb1','FSubpb1');

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
 

%% ============== Subproblem 2: Level Iter ==============
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------ Solving Subpb2 using FEM ------
disp(['--- Start step ',num2str(ALSolveStep),' subproblem 2 ---']); tic;
muvar = 1e-3; udual = 0*FSubpb1; vdual = 0*USubpb1; 
%%JY!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %betavarList = [1e-3,sqrt(1e-5),1e-2,sqrt(1e-3),1e-1,sqrt(1e-1),1e0]*mean(winstepsize/(2^(LevelNo-1))).^2.*mu;
% betavarList = [1e-3,sqrt(1e-5),1e-2,sqrt(1e-3),1e-1,sqrt(1e-1),1e0]*mean(winstepsize).^2.*mu;
% % betavarList = 1e-2 * mean(winstepsize/(2^(LevelNo-1))).^2.*mu;
% alphaList = [ 1e2 ]; alpha2 = 0;
% %Err1 = zeros(length(betavarList),length(alphaList)); Err2 = Err1;
% %%betavar = 0.1*(winstepsize/(2^(LevelNo-1)))^2*mu;  
% %%if LevelNo > 3, betavar = 0.1*(winstepsize/(2^(3-1)))^2*mu; end
% [tempkk,templl] = ndgrid(1:1:length(betavarList),1:1:length(alphaList));
% tempkk = tempkk(:); templl = templl(:); Err1 = 0*tempkk; Err2 = Err1;
% parfor tempm = 1:length(tempkk)
%     tempk = tempkk(tempm); templ = templl(tempm);
% %for tempk = 1:length(betavarList)
% %    for templ = 1:length(alphaList)
%     
%         betavar = betavarList(tempk); alphavar = alphaList(templ)*betavar; alpha2 = 0;
alphavar = 100*betavar; alpha2 = 0;
%         [USubpb2temp,EnrHAndTipEleIndex,EnrTipEleIndex] = Subpb2_adapt(coordinatesFEMLevel2,elementsFEMLevel2,[],neumannLevel2,betavar,mu,...
%             USubpb1,FSubpb1,udual,vdual,alphavar,alpha2,winstepsize, ...
%             CrackOrNot,CrackPath1,CrackPath2,CrackTip,CrackTipOrNot,ALSolveStep); % [] means I don't apply dirichlet/neumann BC here.
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         if CrackOrNot > 0
%             [USubpb2,~] = ComputeCrackDisp(USubpb2temp,x0,y0,coordinatesFEMLevel2,elementsFEMLevel2,EnrTipEleIndex,CrackOrNot,CrackPath1,CrackPath2,CrackTip,CrackTipOrNot);
%             USubpb2 = full(USubpb2); %close all; Plotdisp_show(U0Level2,coordinatesFEMLevel2,elementsFEMLevel2);
%         else, USubpb2 = USubpb2temp; USubpb2 = full(USubpb2); 
%         end
%         [FSubpb2,~,~] = Global_NodalStrainAvg(coordinatesFEMLevel2,elementsFEMLevel2,USubpb2,2,...
%                     winstepsize,LevelNo,EnrHAndTipEleIndex,EnrTipEleIndex,CrackOrNot,CrackPath1,CrackPath2,CrackTip,CrackTipOrNot);
%         % close all; Plotdisp_show(USubpb2,coordinatesFEMLevel2,elementsFEMLevel2);
%         %Plotstrain_show(full(FSubpb2),coordinatesFEMLevel2,elementsFEMLevel2(:,1:4));       
%         %Err1(tempk,templ) = norm(USubpb1-USubpb2,2);
%         %Err2(tempk,templ) = norm(FSubpb1-FSubpb2,2);
%         Err1(tempm) = norm(USubpb1-USubpb2,2); Err2(tempm) = norm(FSubpb1-FSubpb2,2);
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %    end
% end
% Err1 = reshape(Err1,length(betavarList),length(alphaList));   %figure, plot(Err1); surf(Err1);  
% Err2 = reshape(Err2,length(betavarList),length(alphaList));   %figure, plot(Err2); surf(Err2);
% ErrSum = Err1 + Err2*mean(Err2(:))/mean(Err1(:)); figure, plot(ErrSum);  
% [~,indexOfbetavar] = min(ErrSum);
% try
%     [fitobj] = fit(log10(betavarList(indexOfbetavar-1:1:indexOfbetavar+1))',ErrSum(indexOfbetavar-1:1:indexOfbetavar+1),'poly2');
%     p = coeffvalues(fitobj); betavar = 10^(-p(2)/2/p(1));
% catch
%     betavar = betavarList(indexOfbetavar);
% end
 
% ======== Resolve using the new betavar ========
[USubpb2temp,EnrHAndTipEleIndex,EnrTipEleIndex] = Subpb2_adapt(coordinatesFEMLevel2,elementsFEMLevel2,[],neumannLevel2,betavar,muvar,USubpb1,FSubpb1,udual,vdual,alphavar,alpha2, ...
    winstepsize, CrackOrNot,CrackPath1,CrackPath2,CrackTip,CrackTipOrNot,ALSolveStep,ClusterNo); % [] means I don't apply dirichlet/neumann BC here.
ALSub2Time(ALSolveStep) = toc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if CrackOrNot > 0
    [USubpb2,~] = ComputeCrackDisp(USubpb2temp,x0,y0,coordinatesFEMLevel2,elementsFEMLevel2,EnrTipEleIndex,CrackOrNot,CrackPath1,CrackPath2,CrackTip,CrackTipOrNot);
    USubpb2 = full(USubpb2); % close all; Plotdisp_show(USubpb2,coordinatesFEMLevel2,elementsFEMLevel2);
else, USubpb2 = full(USubpb2temp); % close all; Plotdisp_show(USubpb2,coordinatesFEMLevel2,elementsFEMLevel2);
end
% ======== Smooth the displacement field =======
try
  temp = DispFilterSize; temp = DispFilterStd; temp = StrainFilterSize; temp = StrainFilterStd;
catch
 DispFilterSize = 0; DispFilterStd = 0; StrainFilterSize = 0; StrainFilterStd = 0;
end
% USubpb2 = funSmoothDisp_adapt(USubpb2,M,N,x,y);
[FSubpb2,~,~] = Global_NodalStrainAvg(coordinatesFEMLevel2,elementsFEMLevel2,USubpb2,2,...
    winstepsize,LevelNo,EnrHAndTipEleIndex,EnrTipEleIndex,CrackOrNot,CrackPath1,CrackPath2,CrackTip,CrackTipOrNot);

% [FSubpb2,StrainGaussPt,CoordGaussPt] = Global_NodalStrainAvg(coordinatesFEMLevel2,elementsFEMLevel2,USubpb1,2,...
%     winstepsize,LevelNo,EnrHAndTipEleIndex,EnrTipEleIndex,CrackOrNot,CrackPath1,CrackPath2,CrackTip,CrackTipOrNot);
% tritemp = delaunay(CoordGaussPt(:,1),CoordGaussPt(:,2));
% figure,trisurf(tritemp,CoordGaussPt(:,1),CoordGaussPt(:,2),StrainGaussPt(:,1));

% ------- Smooth strain field --------
%FSubpb2 = funSmoothStrain_adapt(FSubpb2,coordinatesFEMLevel2,elementsFEMLevel2,winstepsize,LevelNo,StrainFilterSize,StrainFilterStd);
  

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
        % USubpb2 = funSmoothDisp_adapt(USubpb2,coordinatesFEMLevel2,elementsFEMLevel2,winstepsize,LevelNo,DispFilterSize,DispFilterStd);
        % Plotdisp_show(USubpb2,coordinatesFEMLevel2,elementsFEMLevel2(:,1:4));
        % FSubpb2 = funSmoothStrain_adapt(FSubpb2,coordinatesFEMLevel2,elementsFEMLevel2,winstepsize,LevelNo,StrainFilterSize,StrainFilterStd);
        % Plotstrain_show(FSubpb2,coordinatesFEMLevel2,elementsFEMLevel2(:,1:4));
    end
end
disp(['--- Done! step ',num2str(ALSolveStep),' subproblem 2 ---']);

% ------ save data ------
save(['Subpb2_step',num2str(ALSolveStep)],'USubpb2','FSubpb2'); 
udual = (FSubpb2 - FSubpb1); vdual = (USubpb2 - USubpb1);
save(['uvdual_step',num2str(ALSolveStep)],'udual','vdual');
 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ALSolveStep = 1; UpdateY = 1/tol2; % tol2 = 1e-4;  CrackOrNot = 0; CrackPath1 = [0,0]; CrackPath2 = [0,0]; CrackTip = [0,0]; 
HPar = cell(21,1); for tempj = 1:21, HPar{tempj} = HtempPar(:,tempj); end
% ======== Start ALStep global iteration ======== 
while (ALSolveStep<3)
    ALSolveStep = ALSolveStep + 1;  % Update using the last step
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ============== Subproblem 1 ==============
    [USubpb1, HPar, ALSub1Timetemp, ConvItPerEletemp] = Subpb1_adapt(USubpb2,(0.9*FSubpb1+0.1*FSubpb2),udual,vdual, coordinatesFEMLevel2,...
    fNormalized,gNormalized,Df,winsizeList, CrackOrNot,CrackPath1,CrackPath2,CrackTip, ...
    muvar,betavar,HPar,ALSolveStep,tol,'GaussNewton',winstepsize,ClusterNo);
    ALSub1Time(ALSolveStep) = ALSub1Timetemp; ConvItPerEle(:,ALSolveStep) = ConvItPerEletemp; clear ALSub1Timetemp ConvItPerEletemp;
    FSubpb1 = (0.9*FSubpb1+0.1*FSubpb2); 
    disp(['--- Done! step ',num2str(ALSolveStep),' subproblem 1 ---']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [USubpb1,FSubpb1] = funRemoveOutliers_adapt(coordinatesFEMLevel2,elementsFEMLevel2,USubpb1,FSubpb1,0);
    % ------  Manually find some bad points from Local Subset ICGN step ------
%     figure; show([], elementsFEMLevel2,coordinatesFEMLevel2,USubpb1(1:2:end));
%     title('x-Disp U','fontweight','normal'); set(gca,'fontsize',16); view(2); axis tight; axis equal; colorbar; set(gca,'XTick',[]);
%     figure; show([], elementsFEMLevel2,coordinatesFEMLevel2,USubpb1(2:2:end));
%     title('y-Disp V','fontweight','normal'); set(gca,'fontsize',16); view(2); axis tight; axis equal; colorbar; set(gca,'XTick',[]);
     USubpb1 = funSmoothDisp_adapt(USubpb1,coordinatesFEMLevel2,elementsFEMLevel2,winstepsize,LevelNo,DispFilterSize,DispFilterStd);
        Plotdisp_show(USubpb1,coordinatesFEMLevel2,elementsFEMLevel2(:,1:4));
        FSubpb1 = funSmoothStrain_adapt(FSubpb1,coordinatesFEMLevel2,elementsFEMLevel2,winstepsize,LevelNo,StrainFilterSize,StrainFilterStd);
        Plotstrain_show(FSubpb1,coordinatesFEMLevel2,elementsFEMLevel2(:,1:4));
%     disp('--- Start to manually remove bad points ---')
%     %close all; Plotuv(USubpb1,x0,y0);
%     u = reshape(USubpb1(1:2:end),M,N); v = reshape(USubpb1(2:2:end),M,N);
%     %[u,v,~,Subpb1_BadptRow,Subpb1_BadptCol] = funRemoveOutliers(u,v,[],0,Local_BadptRow,Local_BadptCol);
%     qDICOrNot = 0; Thr0 = 100; [u,v,cc] = funRemoveOutliers(u,v,cc,qDICOrNot,Thr0,Local_BadptRow,Local_BadptCol);
%     disp('--- Remove bad points done ---')
%     USubpb1(1:2:end) = reshape(u,size(elements,1),1); USubpb1(2:2:end) = reshape(v,size(elements,1),1);
%     %close all; Plotuv(USubpb1,x0,y0); 
    % Plotdisp_show(USubpb1-USubpb1_iter1,coordinatesFEMLevel2,elementsFEMLevel2);
% Plotdisp_show(USubpb1-USubpb2,coordinatesFEMLevel2,elementsFEMLevel2);
      save(['Subpb1_step',num2str(ALSolveStep)],'USubpb1','FSubpb1');
     % USubpb1 = funSmoothDisp(USubpb1,coordinatesFEMLevel2,elementsFEMLevel2,[],[],winstepsize,DispFilterStd,DispFilterSize,LevelNo);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ============== Subproblem 2 ==============
    % ------ Solving Subpb2 using FEM ------
    disp(['--- Start step ',num2str(ALSolveStep),' subproblem 2 ---']); tic;
    [USubpb2temp,EnrHAndTipEleIndex,EnrTipEleIndex] = Subpb2_adapt(coordinatesFEMLevel2,elementsFEMLevel2,[],neumannLevel2,betavar,muvar,USubpb1,FSubpb1,udual,vdual,alphavar,alpha2, ...
        winstepsize, CrackOrNot,CrackPath1,CrackPath2,CrackTip,CrackTipOrNot,ALSolveStep,ClusterNo); % [] means I don't apply dirichlet/neumann BC here.
    ALSub2Time(ALSolveStep) = toc;
    disp(['--- Done! step ',num2str(ALSolveStep),' subproblem 2 ---']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    if CrackOrNot > 0
        [USubpb2,~] = ComputeCrackDisp(USubpb2temp,x0,y0,coordinatesFEMLevel2,elementsFEMLevel2,EnrTipEleIndex,CrackOrNot,CrackPath1,CrackPath2,CrackTip,CrackTipOrNot);
        USubpb2 = full(USubpb2); % close all; Plotuv(USubpb2,x0,y0); Plotdisp_show(USubpb2,elementsFEM,coordinatesFEM);
    else
        USubpb2 = USubpb2temp; USubpb2 = full(USubpb2); % close all; Plotuv(USubpb2,x0,y0); 
        % Plotdisp_show(USubpb2,coordinatesFEMLevel2,elementsFEMLevel2);
    end
    
    % ======== Smooth the displacement field =======
    % ULocal = funSmoothDisp(USubpb2,M,N,x,y);
    [FSubpb2,~,~] = Global_NodalStrainAvg(coordinatesFEMLevel2,elementsFEMLevel2,USubpb2,2,...
                    winstepsize,LevelNo,EnrHAndTipEleIndex,EnrTipEleIndex,CrackOrNot,CrackPath1,CrackPath2,CrackTip,CrackTipOrNot);
    
    % ------- Smooth strain field --------
    % FSubpb2 = funSmoothStrain(FSubpb2,coordinatesFEMLevel2,elementsFEMLevel2,winstepsize,StrainFilterSize,StrainFilterStd,LevelNo);

    % ------- Before computing strain, we smooth the displacement field -------
    if LevelNo < 6
        if CrackOrNot > 0
%             % ------- Find elements located on the top/bottom of crack path -------                        
%             [CoordCrackTop,EleCrackTop,CoordCrackBottom,EleCrackBottom] = funCrackTopOrBottom(coordinatesFEMLevel2,elementsFEMLevel2,CrackPath1,CrackPath2,CrackTip);
%             % ------- Before computing strain, we smooth the displacement field -------
%             USubpb2 = funSmoothDispCrack(USubpb2,coordinatesFEMLevel2,elementsFEMLevel2,winstepsize,CoordCrackTop,CoordCrackBottom,LevelNo,DispFilterSize,DispFilterStd);
%             % -------- Plot disp -------- 
%             Plotdisp_show(USubpb2,coordinatesFEMLevel2,elementsFEMLevel2([EleCrackTop;EleCrackBottom],1:4)); 
%             % -------- Compute strain field -------- 
%             FLocal = ComputeCrackStrain(USubpb2,FSubpb2,coordinatesFEMLevel2,elementsFEMLevel2,EleCrackTop,EleCrackBottom,EnrHAndTipEleIndex,EnrTipEleIndex,0,CrackPath1,CrackPath2,CrackTip,0);
%             FSubpb2 = funSmoothStrainCrack(FLocal,coordinatesFEMLevel2,elementsFEMLevel2,winstepsize,CoordCrackTop,CoordCrackBottom,LevelNo,StrainFilterSize,StrainFilterStd);
%             Plotstrain_show(FSubpb2,coordinatesFEMLevel2,elementsFEMLevel2([EleCrackTop;EleCrackBottom],1:4));
        else
            % USubpb2 = funSmoothDisp_adapt(USubpb2,coordinatesFEMLevel2,elementsFEMLevel2,winstepsize,LevelNo,DispFilterSize,DispFilterStd);
            %Plotdisp_show(USubpb2,coordinatesFEMLevel2,elementsFEMLevel2(:,1:4));
            % FSubpb2 = funSmoothStrain_adapt(FSubpb2,coordinatesFEMLevel2,elementsFEMLevel2,winstepsize,LevelNo,StrainFilterSize,StrainFilterStd);
            %Plotstrain_show(FSubpb2,coordinatesFEMLevel2,elementsFEMLevel2(:,1:4));
        end
    end
    % ------ save data ------
    save(['Subpb2_step',num2str(ALSolveStep)],'USubpb2','FSubpb2');
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ==============Update dual variables ===============
    udual = udual + (FSubpb2 - FSubpb1);
    vdual = vdual + (USubpb2 - USubpb1); 
    save(['uvdual_step',num2str(ALSolveStep)],'udual','vdual');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute norm of UpdateY
    USubpb2_Old = load(['Subpb2_step',num2str(ALSolveStep-1)],'USubpb2');
    USubpb2_New = load(['Subpb2_step',num2str(ALSolveStep)],'USubpb2');
    UpdateY = norm((USubpb2_Old.USubpb2 - USubpb2_New.USubpb2), 2)/sqrt(size(USubpb2_Old.USubpb2,1)); UpdateY
    
    USubpb1_Old = load(['Subpb1_step',num2str(ALSolveStep-1)],'USubpb1');
    USubpb1_New = load(['Subpb1_step',num2str(ALSolveStep)],'USubpb1');
    UpdateY2 = norm((USubpb1_Old.USubpb1 - USubpb1_New.USubpb1), 2)/sqrt(size(USubpb1_Old.USubpb1,1)); UpdateY2
    
    if UpdateY < tol2 || UpdateY2 < tol; break; end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UIter = full(USubpb2); FIter = full(FSubpb2); 
udualIter=full(udual); vdualIter=full(vdual); USubpb2tempIter = full(USubpb2temp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------- Plot disp and strain field -------- 
close all; 
if CrackOrNot > 0
    Plotdisp_show(UIter,coordinatesFEMLevel2,elementsFEMLevel2([EleCrackTop;EleCrackBottom],1:4)); 
    Plotstrain_show(FIter,coordinatesFEMLevel2,elementsFEMLevel2([EleCrackTop;EleCrackBottom],1:4));
else
    Plotdisp_show(UIter,coordinatesFEMLevel2,elementsFEMLevel2(:,1:4));
    Plotstrain_show(FIter,coordinatesFEMLevel2,elementsFEMLevel2(:,1:4));
    
          
figure(1); savefig(['fig_L',num2str(LevelNo),'_dispu.fig']);
figure(2); savefig(['fig_L',num2str(LevelNo),'_dispv.fig']);
figure(3); savefig(['fig_L',num2str(LevelNo),'_e11.fig']);
figure(4); savefig(['fig_L',num2str(LevelNo),'_e22.fig']);
figure(5); savefig(['fig_L',num2str(LevelNo),'_e12.fig']);
figure(6); savefig(['fig_L',num2str(LevelNo),'_eshear.fig']);

end
 
% Plotstrain_show(FSubpb1,coordinatesFEMLevel2,elementsFEMLevel2(:,1:4));
disp(['------------ LevelNo ',num2str(LevelNo),' Done ------------ ']);


% ------ Delete temp files ------
for tempi = 1:ALSolveStep
    file_name_Subpb1 = ['Subpb1_step',num2str(tempi),'.mat'];
    file_name_Subpb2 = ['Subpb2_step',num2str(tempi),'.mat'];
    file_name_dual = ['uvdual_step',num2str(tempi),'.mat'];
    delete(file_name_Subpb1); delete(file_name_Subpb2); delete(file_name_dual);
end

betamuIter.betavar = betavar; betamuIter.mu = muvar; betamuIter.alphavar = alphavar;

% ============== Assign values to variables ===============
% ULevel2 = UIter; FLevel2 = FIter; udualLevel2 = udualIter; vdualLevel2 = vdualIter; USubpb2tempLevel2 = USubpb2tempIter;
% ConvItPerEleLevel2 = ConvItPerEle; PassCrackOrNotLevel2 = PassCrackOrNotIter;
% rhoLevel2Vector = rhoIterVector; rhoLevel2Vectortemp1 = rhoIterVectortemp1; rhoLevel2Vectortemp2 = rhoIterVectortemp2;
% U0Level3 = U0IterNext; F0Level3 = F0IterNext; elementsFEMLevel3 = elementsFEMIterNext; coordinatesFEMLevel3 = coordinatesFEMLevel2Next;
% refinedEleIDListLevel2 = refinedEleIDListIter; eleGenerationLevel3 = eleGenerationIterNext;
% dirichletLevel3 = dirichletIterNext(2:end); neumannLevel3 = [];
% ALSolveStepLevel2 = ALSolveStep;
% ALSub1TimeLevel2 = ALSub1Time; ALSub2TimeLevel2 = ALSub2Time; EstimateTimeLevel2 = EstimateTime; MarkTimeLevel2 = MarkTime; RefineTimeLevel2 = RefineTime;




