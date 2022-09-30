%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that do adaptive global DIC 
% in the oder of: (Input: (LevelNo-1)Solve) -> Estimate -> Mark -> Refine -> Solve(Output)
% 
% In this code, variables at (LevelNo-1) are all called "Level1", while
% variables at (LevelNo) are called "Level2".
%
% Author: Jin Yang, Email: jyang526@wisc.edu
% Date: 2019 June.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ULevel2,FLevel2,alphavar,normOfWLevel2,...
    coordinatesFEMLevel2,elementsFEMLevel2,eleGenerationLevel2,irregularEdgeLevel2,dirichletLevel2,neumannLevel2,...
    thetaDorfler,TimeICGNLevel2,EstimateTime,MarkTime,RefineTime,...
    rhoLevel1Vector,rhoLevel1Vectortemp1,rhoLevel1Vectortemp2,EnrHAndTipEleIndex,EnrTipEleIndex] ...
    = funGBMSIterSq(ULevel1,FLevel1,coordinatesFEM,elementsFEM,winstepsize,ClusterNo,...
    eleGeneration,irregularEdge,dirichlet,neumann,LevelNo,fNormalized,gNormalized,winsize,...
    CrackOrNot,CrackPath1,CrackPath2,CrackTip,CrackTipOrNot,EnrHAndTipEleIndex,EnrTipEleIndex,...
    DispFilterSize,DispFilterStd,StrainFilterSize,StrainFilterStd,tol,alphavarLevel1)

%elementsFEM = elementsFEMLevel1;  coordinatesFEM = coordinatesFEMLevel1;
%neumann = neumannLevel1; dirichlet = dirichletLevel1;
alphavar = alphavarLevel1; GaussPtOrder = 2;
imgPyramidUnit = 1; Df = funImgGradient(fNormalized,gNormalized); 

%% ============== Estimate: Level1 ==============
% ------ Compute rho to find where to refine ------
rhoLevel1Vector = zeros(size(elementsFEM,1),1); rhoLevel1Vectortemp1 = rhoLevel1Vector; rhoLevel1Vectortemp2 = rhoLevel1Vector;
ULevel1temp = [ULevel1;0;0]; imgPyramidUnit = 1; 
try
    hbar = parfor_progressbar(size(elementsFEM,1),'Please wait for a posterior error estimate!'); tic
    parfor j = 1:size(elementsFEM,1)
    % ---- Find the element nodal indices ----
    nodalindextemp = [2*elementsFEM(j,1)-1 2*elementsFEM(j,1) 2*elementsFEM(j,2)-1 2*elementsFEM(j,2) ...
            2*elementsFEM(j,3)-1 2*elementsFEM(j,3) 2*elementsFEM(j,4)-1 2*elementsFEM(j,4) (size(ULevel1,1)+1)*ones(1,8)];
    for tempi = 5:8    
        if elementsFEM(j,tempi) ~= 0
            nodalindextemp(2*tempi-1) = 2*elementsFEM(j,tempi)-1; nodalindextemp(2*tempi) = 2*elementsFEM(j,tempi);
        end
    end 
    % ---- Find element corner points ----
    point1x = coordinatesFEM(elementsFEM(j,1),1); point3x = coordinatesFEM(elementsFEM(j,3),1); lengthOfElement = (point3x-point1x)/winstepsize;
    % ---- Compute error estimator ----
    rhoLevel1Vectortemp1(j) = (lengthOfElement)^2 * aPostErrEstIntSq(coordinatesFEM,elementsFEM,j,fNormalized,gNormalized,Df,ULevel1temp(nodalindextemp));
    rhoLevel1Vectortemp2(j) = (lengthOfElement) * alphavar * aPostErrEstJumpSq(coordinatesFEM,elementsFEM,j,eleGeneration,winstepsize,ULevel1);
    rhoLevel1Vector(j) = rhoLevel1Vectortemp1(j) + rhoLevel1Vectortemp2(j);
    hbar.iterate(1);
    end 
catch
    hbar = waitbar(0,'Please wait for a posterior error estimate!');
    for j = 1:size(elementsFEM,1)
    % ---- Find the element nodal indices ----
    nodalindextemp = [2*elementsFEM(j,1)-1 2*elementsFEM(j,1) 2*elementsFEM(j,2)-1 2*elementsFEM(j,2) ...
            2*elementsFEM(j,3)-1 2*elementsFEM(j,3) 2*elementsFEM(j,4)-1 2*elementsFEM(j,4) (size(ULevel1,1)+1)*ones(1,8)];
    for tempi = 5:8    
        if elementsFEM(j,tempi) ~= 0
            nodalindextemp(2*tempi-1) = 2*elementsFEM(j,tempi)-1; nodalindextemp(2*tempi) = 2*elementsFEM(j,tempi);
        end
    end 
    % ---- Find element corner points ----
    point1x = coordinatesFEM(elementsFEM(j,1),1); point3x = coordinatesFEM(elementsFEM(j,3),1); lengthOfElement = (point3x-point1x)/winstepsize;
    % ---- Compute error estimator ----
    rhoLevel1Vectortemp1(j) = (lengthOfElement)^2 * aPostErrEstIntSq(coordinatesFEM,elementsFEM,j,fNormalized,gNormalized,Df,ULevel1temp(nodalindextemp));
    rhoLevel1Vectortemp2(j) = (lengthOfElement) * alphavar * aPostErrEstJumpSq(coordinatesFEM,elementsFEM,j,eleGeneration,winstepsize,ULevel1);
    rhoLevel1Vector(j) = rhoLevel1Vectortemp1(j) + rhoLevel1Vectortemp2(j);
    waitbar(j/size(elementsFEM,1));
    end 
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

% ======= Compute SSD Error and Penalty Error =======
[ErrSSD] = funSSDErr(coordinatesFEM,elementsFEM,fNormalized,gNormalized,ULevel1,0,0);
sum(ErrSSD(:))
% % ------- Plot SSD ---------
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
    thetaDorflerList(LevelNo) = thetaDorfler; rhoLevel1SumM = 0; rhoLevel1SumAll = sum(rhoLevel1Vector); rhoLevel1SumMIndex = 0;
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
    figure; patch(Sqx,Sqy,Sqc,'EdgeColor','none'); title('Marked ele (yellow)','fontweight','normal'); set(gca,'fontsize',16); axis tight; axis equal;  

     savefig(['fig_L',num2str(LevelNo+1),'_mark.fig']); % JY!!!
    
    
   %% ============== Refine: Level1 ==============
    if thetaDorfler>0 && size(coordinatesFEM,1)>1e4  
    
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
        
     
    end
    
    

    % ====== Show refined mesh ======
    Sqx = zeros(4,size(elementsFEMLevel2,1)); Sqy = zeros(4,size(elementsFEMLevel2,1)); Sqc = zeros(1,size(elementsFEMLevel2,1));
    for j = 1:size(elementsFEMLevel2,1)
        Sqx(1:4,j) = coordinatesFEMLevel2(elementsFEMLevel2(j,1:4),1);
        Sqy(1:4,j) = coordinatesFEMLevel2(elementsFEMLevel2(j,1:4),2); Sqc(j) =  1;
    end
    figure; patch(Sqx,Sqy,Sqc);  axis equal; axis off; title('Level2 Mesh','fontweight','normal'); set(gca,'fontsize',16);

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
LevelNo = LevelNo + 1;
% 
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%% ============== Level2 ==============        
PassCrackOrNotIter = zeros(size(elementsFEMLevel2,1),1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ULevel2, normOfWLevel2, TimeICGNLevel2] = funGlobalICGN_adapt(coordinatesFEMLevel2,elementsFEMLevel2,Df,fNormalized,gNormalized,U0Level2,alphavar,tol,ClusterNo);
[FLevel2] = funGlobal_NodalStrainAvg(coordinatesFEMLevel2,elementsFEMLevel2,ULevel2,GaussPtOrder);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[USubpb1,FSubpb1] = funRemoveOutliers_adapt(coordinatesFEMLevel2,elementsFEMLevel2,USubpb1,FSubpb1,0);
% close all; Plotdisp_show(USubpb1,coordinatesFEMLevel2,elementsFEMLevel2);
% Plotstrain_show(FSubpb1,coordinatesFEMLevel2,elementsFEMLevel2);
% save(['Subpb1_step',num2str(ALSolveStep)],'USubpb1','FSubpb1');

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
% [FSubpb2,~,~] = Global_NodalStrainAvg(coordinatesFEMLevel2,elementsFEMLevel2,USubpb2,2,...
%     winstepsize,LevelNo,EnrHAndTipEleIndex,EnrTipEleIndex,CrackOrNot,CrackPath1,CrackPath2,CrackTip,CrackTipOrNot);

% [FSubpb2,StrainGaussPt,CoordGaussPt] = Global_NodalStrainAvg(coordinatesFEMLevel2,elementsFEMLevel2,USubpb1,2,...
%     winstepsize,LevelNo,EnrHAndTipEleIndex,EnrTipEleIndex,CrackOrNot,CrackPath1,CrackPath2,CrackTip,CrackTipOrNot);
% tritemp = delaunay(CoordGaussPt(:,1),CoordGaussPt(:,2));
% figure,trisurf(tritemp,CoordGaussPt(:,1),CoordGaussPt(:,2),StrainGaussPt(:,1));

% ------- Smooth strain field --------
%FLevel2 = funSmoothStrain_adapt(FLevel2,coordinatesFEMLevel2,elementsFEMLevel2,winstepsize,LevelNo,StrainFilterSize,StrainFilterStd);
  

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
          Plotdisp_show(ULevel2,coordinatesFEMLevel2,elementsFEMLevel2(:,1:4));
        % FSubpb2 = funSmoothStrain_adapt(FSubpb2,coordinatesFEMLevel2,elementsFEMLevel2,winstepsize,LevelNo,StrainFilterSize,StrainFilterStd);
          Plotstrain_show(FLevel2,coordinatesFEMLevel2,elementsFEMLevel2(:,1:4));
        
          figure(1); savefig(['fig_L',num2str(LevelNo),'_dispu.fig']);
          figure(2); savefig(['fig_L',num2str(LevelNo),'_dispv.fig']);
          figure(3); savefig(['fig_L',num2str(LevelNo),'_e11.fig']);
          figure(4); savefig(['fig_L',num2str(LevelNo),'_e22.fig']);
          figure(5); savefig(['fig_L',num2str(LevelNo),'_e12.fig']);
          figure(6); savefig(['fig_L',num2str(LevelNo),'_eshear.fig']);
 
    end
end
disp(['--- Done! step ',num2str(LevelNo),' ---']);
 