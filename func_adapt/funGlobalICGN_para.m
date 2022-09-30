%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function FE-based Global DVC ICGN-code                   %
% Object: to find deformation field using global methods   %
% Author: Jin Yang                                         %
% Last date modified: 2019.03                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [U,normOfW,TimeICGN] = funGlobalICGN_para(coordinatesFEM,elementsFEM,Df,Img1,Img2,U,alpha,tol,clusterNo)
    
DIM = 2; NodesPerEle = 4; % Using cubic elements
FEMSize = DIM*size(coordinatesFEM,1); % FE-system size
winsize = (coordinatesFEM(2,1)-coordinatesFEM(1,1))*ones(1,DIM);

DfDx = Df.DfDx; DfDy = Df.DfDy;  
DfAxis = Df.DfAxis; DfDxStartx = DfAxis(1); DfDxStarty = DfAxis(3);  

for stepwithinwhile = 1:100
    
   tic; disp(['--- Global IC-GN iteration step',num2str(stepwithinwhile),' ---']);
   
% ====== Initialize stiffness matrix at first iteration ======
%    if (stepwithinwhile==1)
%        INDEXAI = []; INDEXAJ = []; INDEXAVAL = []; INDEXAREG = []; 
%        %A = sparse(FEMSize,FEMSize);
%    end
%    INDEXBI = []; INDEXBVAL = []; %clear b; b = sparse(FEMSize,1);
      
   % ------ Define ksi, eta and zeta list ------
   if stepwithinwhile==1
   ksiList = -1:2/winsize(1):1; etaList = -1:2/winsize(2):1;  
   [ksiMat,etaMat] = ndgrid(ksiList,etaList);
   
   NMat = cell(4,1);
   NMat{1} = 1/4*(1-ksiMat).*(1-etaMat); NMat{2} = 1/4*(1+ksiMat).*(1-etaMat);
   NMat{3} = 1/4*(1+ksiMat).*(1+etaMat); NMat{4} = 1/4*(1-ksiMat).*(1+etaMat);
   end
   
   if clusterNo== 0 || clusterNo== 1
       hbar = waitbar(0, ['Global ICGN iteartion step: ',num2str(stepwithinwhile)]);
   else
       hbar=parfor_progressbar(size(elementsFEM,1),['Global ICGN iteartion step: ',num2str(stepwithinwhile)]);
   end
   
   % ============= Each element, assemble stiffness matrix ============
   parfor j = 1: size(elementsFEM,1) % j is the element index
        
        if clusterNo== 0 || clusterNo== 1
            waitbar(j/size(elementsFEM,1));
        else
            hbar.iterate(1);
        end
        
        tempA = zeros(DIM*NodesPerEle,DIM*NodesPerEle); tempb = tempA(:,1);
        
        % ------ Find corner pts ------
        pt1xyz = coordinatesFEM(elementsFEM(j,1),:); pt2xyz = coordinatesFEM(elementsFEM(j,2),:);
        pt3xyz = coordinatesFEM(elementsFEM(j,3),:); pt4xyz = coordinatesFEM(elementsFEM(j,4),:);
        
        pt1x = pt1xyz(1); pt1y = pt1xyz(2);  pt2x = pt2xyz(1); pt2y = pt2xyz(2);  
        pt3x = pt3xyz(1); pt3y = pt3xyz(2);  pt4x = pt4xyz(1); pt4y = pt4xyz(2);  
         
        % ------ Calculate ksi and eta ------
        % lMatrix = [point1x*point1y point1x point1y 1;
        %     point2x*point2y point2x point2y 1;
        %     point3x*point3y point3x point3y 1;
        %     point4x*point4y point4x point4y 1];
        
        % ------ Find linear interpolation coefficients ------
        % lb = [-1;1;1;-1]; lCoeff = linsolve(lMatrix,lb);
        % mb = [-1;-1;1;1]; mCoeff = linsolve(lMatrix,mb);
        
        % ------ Find the element nodal indices ------
        tp = ones(1,DIM);
        tempIndexU = 2*elementsFEM(j,[tp,2*tp,3*tp,4*tp]);
        tempIndexU(1:2:end) = tempIndexU(1:2:end)-1;  % size of tempIndexU: 1*8
        
        [ptOfxAll,ptOfyAll] = ndgrid(pt1x:pt3x,pt1y:pt3y);
          
        % U1Mat = U(tempIndexU(3*1-2))*ones(winsize+ones(1,3)); U2Mat = U(tempIndexU(3*2-2))*ones(winsize+ones(1,3));
        % U3Mat = U(tempIndexU(3*3-2))*ones(winsize+ones(1,3)); U4Mat = U(tempIndexU(3*4-2))*ones(winsize+ones(1,3));
        % U5Mat = U(tempIndexU(3*5-2))*ones(winsize+ones(1,3)); U6Mat = U(tempIndexU(3*6-2))*ones(winsize+ones(1,3));
        % U7Mat = U(tempIndexU(3*7-2))*ones(winsize+ones(1,3)); U8Mat = U(tempIndexU(3*8-2))*ones(winsize+ones(1,3));
        tempUMat = zeros(winsize+ones(1,2)); tempVMat = tempUMat;  
        for tempk = 1:NodesPerEle
            tempUMat = tempUMat + (U(tempIndexU(2*tempk-1))*ones(winsize+ones(1,2))).*NMat{tempk};
            tempVMat = tempVMat + (U(tempIndexU(2*tempk-0))*ones(winsize+ones(1,2))).*NMat{tempk};
        end
        %ptOfxMat = ptOfxAll + tempUMat; ptOfyMat = ptOfyAll + tempVMat; ptOfzMat = ptOfzAll + tempWMat;
        
        tempg = ba_interp2(Img2, ptOfyAll+tempVMat, ptOfxAll+tempUMat,'cubic');
         
        ptOfxAll = ptOfxAll(:); ptOfyAll = ptOfyAll(:); 
        %for ptOfx = pt1x:pt7x
        %    for ptOfy = pt1y:pt7y
        %        for ptOfz = pt1z:pt7z
        for tempjj = 1:length(ptOfxAll)
            
            %waitbar(tempjj/length(ptOfxAll));
            
            ptOfx = ptOfxAll(tempjj); ptOfy = ptOfyAll(tempjj);  
             
                    % ------ Calculate ksi, eta and zeta ------
                    ksi = ksiMat(tempjj); eta = etaMat(tempjj);  
                    %ksi = [ptOfx*ptOfy*ptOfz, ptOfx*ptOfy, ptOfy*ptOfz, ptOfz*ptOfx, ptOfx, ptOfy, ptOfz, 1]*lCoeff;
                    %eta = [ptOfx*ptOfy*ptOfz, ptOfx*ptOfy, ptOfy*ptOfz, ptOfz*ptOfx, ptOfx, ptOfy, ptOfz, 1]*mCoeff;
                    %zeta = [ptOfx*ptOfy*ptOfz, ptOfx*ptOfy, ptOfy*ptOfz, ptOfz*ptOfx, ptOfx, ptOfy, ptOfz, 1]*nCoeff;
                    
                    % ------ Calculate N matrix ------
                    N1 = NMat{1}(tempjj); N2 = NMat{2}(tempjj); N3 = NMat{3}(tempjj); N4 = NMat{4}(tempjj);
                    % N1 = (1-ksi)*(1-eta)*0.25;
                    % N2 = (1+ksi)*(1-eta)*0.25;
                    % N3 = (1+ksi)*(1+eta)*0.25;
                    % N4 = (1-ksi)*(1+eta)*0.25;
                    
                    % ------ Generate [N] shape function matrix ------
                    tpN1 = diag([N1,N1 ]); tpN2 = diag([N2,N2 ]); tpN3 = diag([N3,N3 ]); tpN4 = diag([N4,N4 ]);
                     
                    NOrig = [tpN1,tpN2,tpN3,tpN4]; N = NOrig;
                    NDiag = diag([N1*ones(1,DIM),N2*ones(1,DIM),N3*ones(1,DIM),N4*ones(1,DIM)]);
                    
                    % ------ Build J matrix ------
                    % Comment: I didn't change Jacobian matrix J when enriched
                    % functions are added.
                    J = [funDN1Dksi(ksi,eta),funDN2Dksi(ksi,eta),funDN3Dksi(ksi,eta),funDN4Dksi(ksi,eta);
                         funDN1Deta(ksi,eta),funDN2Deta(ksi,eta),funDN3Deta(ksi,eta),funDN4Deta(ksi,eta)] * ...
                        [pt1x,pt1y;pt2x,pt2y;pt3x,pt3y;pt4x,pt4y];
                    
                    Jacobian = det(J);
                    InvJ = inv(J); 
                    
                    % ------ Compute [DN] matrix ------
                    DN = [InvJ zeros(2,2); zeros(2,2) InvJ] * ...
                    [funDN1Dksi(ksi,eta) 0 funDN2Dksi(ksi,eta) 0 funDN3Dksi(ksi,eta) 0 funDN4Dksi(ksi,eta) 0;
                     funDN1Deta(ksi,eta) 0 funDN2Deta(ksi,eta) 0 funDN3Deta(ksi,eta) 0 funDN4Deta(ksi,eta) 0;
                     0 funDN1Dksi(ksi,eta) 0 funDN2Dksi(ksi,eta) 0 funDN3Dksi(ksi,eta) 0 funDN4Dksi(ksi,eta);
                     0 funDN1Deta(ksi,eta) 0 funDN2Deta(ksi,eta) 0 funDN3Deta(ksi,eta) 0 funDN4Deta(ksi,eta)];
                  
                    % ------ Here approximate Dg(x+u)=Df(x) ------
                    DfEle = [DfDx(ptOfx-DfDxStartx, ptOfy-DfDxStarty); 
                            DfDy(ptOfx-DfDxStartx, ptOfy-DfDxStarty)];
                  
                    % ------ Only assemble stiffness in the first step ------
                    if (stepwithinwhile==1)
                        %A(tempIndexU,tempIndexU) = A(tempIndexU,tempIndexU) + (N'*DfEle)*(N'*DfEle)' + alpha*(DN')*DN; 
                        tempA = tempA + (N'*DfEle)*(N'*DfEle)' + alpha*(DN')*DN;
                    end
                    
                    % ------ Construct b vector ------
                    %temp1 = [ptOfx;ptOfy;ptOfz] + N*U(tempIndexU);
                    %temp2 = ((Img1(ptOfx,ptOfy,ptOfz) - fungInterpolation_g3(temp1(1),temp1(2),temp1(3),Img2(floor(temp1(1))-1:floor(temp1(1))+2, ...
                    %        floor(temp1(2))-1:floor(temp1(2))+2, floor(temp1(3))-1:floor(temp1(3))+2))) * (N'*DfEle));
                    temp2 = ((Img1(ptOfx,ptOfy) - tempg(tempjj)) * (N'*DfEle));
                    %b(tempIndexU) = b(tempIndexU) + tempb - (alpha*(DN')*DN)*U(tempIndexU);
                    tempb = tempb + temp2 - (alpha*(DN')*DN)*U(tempIndexU); 
               %end
           %end
        %end
        end
   
        % Store tempA and tempb
        if (stepwithinwhile==1) 
            %[IndexAXX,IndexAYY] = ndgrid(tempIndexU,tempIndexU); 
            %INDEXAI = [INDEXAI;IndexAXX(:)]; INDEXAJ = [INDEXAJ;IndexAYY(:)]; INDEXAVAL = [INDEXAVAL;tempA(:)]; %INDEXAREG = [INDEXAREG;tempAreg(:)];
            [IndexAXX,IndexAYY] = ndgrid(tempIndexU,tempIndexU);
            INDEXAIpar{j}=IndexAXX(:); INDEXAJpar{j}=IndexAYY(:); INDEXAVALpar{j}=tempA(:);
        end
        %INDEXBI = [INDEXBI;tempIndexU(:)]; INDEXBVAL = [INDEXBVAL;tempb(:)];
        INDEXBIpar{j}=tempIndexU(:); INDEXBVALpar{j}=tempb(:);   
   
   end
   
   close(hbar);
   
%    if (stepwithinwhile==1)
%         % A = sparse(DIM*FEMSize+NodesPerEle*DIM, DIM*FEMSize+NodesPerEle*DIM);
%         % b = sparse(DIM*FEMSize+NodesPerEle*DIM,1);
%         A = sparse(INDEXAI,INDEXAJ,INDEXAVAL,FEMSize,FEMSize);
%         %Areg = sparse( INDEXAI, INDEXAJ, INDEXAREG , DIM*FEMSize+NodesPerEle*DIM, DIM*FEMSize+NodesPerEle*DIM);
%     end
%     b = sparse(INDEXBI,ones(length(INDEXBI),1),INDEXBVAL,FEMSize,1);
 
     A = sparse(FEMSize, FEMSize); b = sparse(FEMSize,1);
    for eleInd = 1:size(elementsFEM,1)
        A = A + sparse(INDEXAIpar{eleInd}, INDEXAJpar{eleInd}, INDEXAVALpar{eleInd}, FEMSize,FEMSize) ;
        b = b + sparse(INDEXBIpar{eleInd},ones(length(INDEXBIpar{eleInd}),1),INDEXBVALpar{eleInd}, FEMSize,1) ;
    end
     
%     coordsIndexInvolved = unique([0;elementsFEM(:)]); % Need modification for triangle elementsFEM
%     
%     % In adaptive mesh, using the following code:
%     UIndexInvolved = zeros(DIM*(length(coordsIndexInvolved)-1),1);
%     
%     % Not including the first 0-th entry
%     for tempi = 1:(size(coordsIndexInvolved,1)-1)
%         UIndexInvolved(3*tempi-2:3*tempi) = [3*coordsIndexInvolved(tempi+1)-2; ...
%                     3*coordsIndexInvolved(tempi+1)-1; 3*coordsIndexInvolved(tempi+1)];
%     end
%      
%     % ========= Set Dirichlet and Neumann boundary conditions =========
%     if isempty(dirichlet) ~= 1
%         dirichlettemp = [3*dirichlet(:); 3*dirichlet(:)-1; 3*dirichlet(:)-2];
%     else
%         dirichlettemp = [];
%     end
%     if isempty(neumann) ~= 1
%         neumanntemp = [3*neumann(:,1); 3*neumann(:,1)-1; 3*neumann(:,1)-2; 3*neumann(:,2); 3*neumann(:,2)-1; 3*neumann(:,2)-2];
%     else
%         neumanntemp = [];
%     end
%     FreeNodes = setdiff(UIndexInvolved,unique([dirichlettemp]));
    
    % ========= Neumann conditions ===========
    % Last step boundary condition force
    % BCForce = -Global_NodalStrainAvg(coordinatesFEM,elementsFEM,Uhat);
    % for tempj = 1:size(neumann,1)
    %     b(2*neumann(tempj,1:2)-1) = b(2*neumann(tempj,1:2)-1) + norm(coordinatesFEM(neumann(tempj,1),:)-coordinatesFEM(neumann(tempj,2),:)) ...
    %         * ( BCForce(4*neumann(tempj,1:2)-3) * neumann(tempj,3) + BCForce(4*neumann(tempj,1:2)-1) * neumann(tempj,4) );
    %     b(2*neumann(tempj,1:2))   = b(2*neumann(tempj,1:2)) + norm(coordinatesFEM(neumann(tempj,1),:)-coordinatesFEM(neumann(tempj,2),:)) ...
    %         * ( BCForce(4*neumann(tempj,1:2)-1) * neumann(tempj,3) + BCForce(4*neumann(tempj,1:2)) * neumann(tempj,4) );
    % end
    
    % ========= Dirichlet conditions ==========
%     UhatOld = Uhat;
%     UhatNew = sparse(DIM*FEMSize + NodesPerEle*DIM, 1);
%     UhatNew(3*unique(dirichlet)) = U(3*unique(dirichlet));
%     UhatNew(3*unique(dirichlet)-1) = U(3*unique(dirichlet)-1);
%     UhatNew(3*unique(dirichlet)-2) = U(3*unique(dirichlet)-2);
     
    % ========= Solve FEM problem ===========
    % W = A\b;
    W = cgsolve(A,b,tol);
    
    normW = norm(W)/sqrt(size(W,1))
    normOfW(stepwithinwhile) = normW;
    TimeICGN(stepwithinwhile) = toc;
    toc
    U = reshape(U,length(U),1); W = reshape(W,length(W),1);
    
    if stepwithinwhile == 1
        normWOld = normW*10;
    else
        normWOld = normOfW(stepwithinwhile-1);
    end
    
    if (normW < tol) || ((normW/normWOld > 1-tol) && (normW/normWOld < 1))
        U = U + W;
        break;
    elseif (normW >= tol  && normW < (1/tol)) % || ((normW/normWOld >= 1) && (normW/normWOld < 100)))
        U = U + W;
    else
        warning('Get diverged in Global_ICGN!!!')
        break;
    end
    
     
    
end

end



%% ========= subroutines for  FEM shape function derivatives ========
function DN1x=funDN1Dksi(ksi,eta)
DN1x = -(1-eta)/4 ;
end
function DN1y=funDN1Deta(ksi,eta)
DN1y =  -(1-ksi)/4 ;
end
function DN2x=funDN2Dksi(ksi,eta)
DN2x =  (1-eta)/4 ;
end
function DN2y=funDN2Deta(ksi,eta)
DN2y =  -(1+ksi)/4 ;
end
function DN3x=funDN3Dksi(ksi,eta)
DN3x = (1+eta)/4 ;
end
function DN3y=funDN3Deta(ksi,eta)
DN3y =  (1+ksi)/4 ;
end
function DN4x=funDN4Dksi(ksi,eta)
DN4x = -(1+eta)/4 ;
end
function DN4y=funDN4Deta(ksi,eta)
DN4y = (1-ksi)/4 ;
end
