%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function adaptive mesh FE-based Global DVC ICGN-code     %
% Object: to find deformation field using global methods   %
% Author: Jin Yang                                         %
% Last date modified: 2019.03                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [U,normOfW,TimeICGN] = funGlobalICGN_adapt_para(coordinatesFEM,elementsFEM,Df,f,g,U,alpha,tol,clusterNo)
  
DIM = 2; NodesPerEle = 8; % Using cubic elements
FEMSize = DIM*2*size(coordinatesFEM,1)+DIM; % FE-system size
winsize = (coordinatesFEM(2,1)-coordinatesFEM(1,1))*ones(1,DIM);

DfDx = Df.DfDx; DfDy = Df.DfDy;  
DfAxis = Df.DfAxis; DfDxStartx = DfAxis(1); DfDxStarty = DfAxis(3);  

U = [U;zeros(DIM,1)]; dirichlet = []; neumann = [];

for stepwithinwhile = 1:100
    
    tic; disp(['--- Global IC-GN iteration step',num2str(stepwithinwhile),' ---']);
    
%     if (stepwithinwhile==1)
%         INDEXAI = []; INDEXAJ = []; INDEXAVAL = []; INDEXAREG = [];
%         %A = sparse(2*size(coordinatesFEM,1)+2,2*size(coordinatesFEM,1)+2);
%         % A = sparse(2*size(unique(elementsFEMFEMLevel4_MS),1), 2*size(unique(elementsFEMFEMLevel4_MS),1));
%     end
%     INDEXBI = []; INDEXBVAL = []; %clear b;b = sparse(2*size(coordinatesFEM,1)+2,1);
      
    % ------ Define ksi, eta and zeta list ------
    if stepwithinwhile==1
        ksiList = -1:2/winsize(1):1; etaList = -1:2/winsize(2):1;
        [ksiMat,etaMat] = ndgrid(ksiList,etaList);
        
    end

    % =====================================================================
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
        
        if elementsFEM(j,1)~= 0 % Only used in the old version codes, not in the new codes
            
            % ------ Find corner pts ------
            pt1xyz = coordinatesFEM(elementsFEM(j,1),:); pt2xyz = coordinatesFEM(elementsFEM(j,2),:);
            pt3xyz = coordinatesFEM(elementsFEM(j,3),:); pt4xyz = coordinatesFEM(elementsFEM(j,4),:);
            
            pt1x = pt1xyz(1); pt1y = pt1xyz(2);  pt2x = pt2xyz(1); pt2y = pt2xyz(2);
            pt3x = pt3xyz(1); pt3y = pt3xyz(2);  pt4x = pt4xyz(1); pt4y = pt4xyz(2);
            
            % ------ Find mid points 5/6/7/8 -------
            if elementsFEM(j,5) ~= 0
                pt5x = coordinatesFEM(elementsFEM(j,5),1); pt5y = coordinatesFEM(elementsFEM(j,5),2);
            else, pt5x = 0; pt5y = 0; end
            if elementsFEM(j,6) ~= 0
                pt6x = coordinatesFEM(elementsFEM(j,6),1); pt6y = coordinatesFEM(elementsFEM(j,6),2);
            else, pt6x = 0; pt6y = 0; end
            if elementsFEM(j,7) ~= 0
                pt7x = coordinatesFEM(elementsFEM(j,7),1); pt7y = coordinatesFEM(elementsFEM(j,7),2);
            else, pt7x = 0; pt7y = 0; end
            if elementsFEM(j,8) ~= 0
                pt8x = coordinatesFEM(elementsFEM(j,8),1); pt8y = coordinatesFEM(elementsFEM(j,8),2);
            else, pt8x = 0; pt8y = 0; end
             
            % ------ Calculate ksi and eta --------
            % lMatrix = [pt1x*pt1y pt1x pt1y 1; 
            %            pt2x*pt2y pt2x pt2y 1; 
            %            pt3x*pt3y pt3x pt3y 1; 
            %            pt4x*pt4y pt4x pt4y 1];
                    
            % ------ Find the element nodal indices ------
            % lb = [-1;1;1;-1]; l = linsolve(lMatrix,lb);
            % mb = [-1;-1;1;1]; m = linsolve(lMatrix,mb);
             
            % =====================================================================
            % ------- Navier-Lame elasticity regularization --------
            % MatrixGrad = [-1 0 0 0 0 0 0 1; 0 1 0 0 0 0 0 -1; 0 -1 0 0 1 0 0 0; 0 0 1 0 -1 0 0 0;
            %               0 0 1 0 0 -1 0 0; 0 0 0 -1 0 1 0 0; 0 0 0 1 0 0 -1 0; -1 0 0 0 0 0 1 0];
            % MatrixGradUpdate = MatrixGrad(:,1:8);
            % 
            % if elementsFEM(j,5) == 0
            %     MatrixGradUpdate(:,2) = MatrixGradUpdate(:,2) + 0.5*MatrixGrad(:,5);
            %     MatrixGradUpdate(:,3) = MatrixGradUpdate(:,3) + 0.5*MatrixGrad(:,5);
            %     MatrixGradUpdate(:,5) = 0*MatrixGradUpdate(:,5);
            % end
            % if elementsFEM(j,6) == 0
            %     MatrixGradUpdate(:,3) = MatrixGradUpdate(:,3) + 0.5*MatrixGrad(:,6);
            %     MatrixGradUpdate(:,4) = MatrixGradUpdate(:,4) + 0.5*MatrixGrad(:,6);
            %     MatrixGradUpdate(:,6) = 0*MatrixGradUpdate(:,6);
            % end
            % if elementsFEM(j,7) == 0
            %     MatrixGradUpdate(:,4) = MatrixGradUpdate(:,4) + 0.5*MatrixGrad(:,7);
            %     MatrixGradUpdate(:,1) = MatrixGradUpdate(:,1) + 0.5*MatrixGrad(:,7);
            %     MatrixGradUpdate(:,7) = 0*MatrixGradUpdate(:,7);
            % end
            % if elementsFEM(j,8) == 0
            %     MatrixGradUpdate(:,1) = MatrixGradUpdate(:,1) + 0.5*MatrixGrad(:,8);
            %     MatrixGradUpdate(:,2) = MatrixGradUpdate(:,2) + 0.5*MatrixGrad(:,8);
            %     MatrixGradUpdate(:,8) = 0*MatrixGradUpdate(:,8);
            % end
            % 
            % % [row,col] = find(elementsFEM(j,5:8)~=(size(coordinatesFEM,1)+1));
            % % for tempi = 1:length(col)
            % %     MatrixGradUpdate(:,4+tempi) = MatrixGrad(:,4+col(tempi));
            % % end
            % 
            % % MatrixGradUpdate(:,1) = MatrixGradUpdate(:,1) + 0.5*MatrixGrad(:,7) + 0.5*MatrixGrad(:,8);
            % % MatrixGradUpdate(:,2) = MatrixGradUpdate(:,2) + 0.5*MatrixGrad(:,8) + 0.5*MatrixGrad(:,5);
            % % MatrixGradUpdate(:,3) = MatrixGradUpdate(:,3) + 0.5*MatrixGrad(:,5) + 0.5*MatrixGrad(:,6);
            % % MatrixGradUpdate(:,4) = MatrixGradUpdate(:,4) + 0.5*MatrixGrad(:,6) + 0.5*MatrixGrad(:,7);
            % 
            % Matrix1 = MatrixGradUpdate'*MatrixGradUpdate;
            % Matrix1Update = zeros(2*size(Matrix1,1),2*size(Matrix1,2));
            % Matrix1Update(1:2:end,1:2:end) = Matrix1;
            % Matrix1Update(2:2:end,2:2:end) = Matrix1;
            % 
            % % Matrix2 = zeros(2*size(MatrixGradUpdate,1),2*size(MatrixGradUpdate,2));
            % % Matrix2(1:2:end,1:2:end) = MatrixGradUpdate;
            % % Matrix2(2:2:end,2:2:end) = MatrixGradUpdate;
            % % Matrix2Update = 0.25*Matrix2'*diag([1 0 1 0 0 1 0 1 1 0 1 0 0 1 0 1])*Matrix2;
            % 
            % % Lame elasticity constants
            % mu = alpha*1; lamda = mu;
            % =====================================================================
            
            % ------ Find the element nodal indices ------
            tp = ones(1,DIM);
            tempIndexU = 2*elementsFEM(j,[tp,2*tp,3*tp,4*tp,5*tp,6*tp,7*tp,8*tp]);
            tempIndexU(1:2:end) = tempIndexU(1:2:end)-1;  % size of tempIndexU: 1*16
            % We don't want temp <= 0, instead, put them to the end
            for tempi = 5:8
                if elementsFEM(j,tempi) == 0
                    tempIndexU(2*tempi-1) = 2*(size(coordinatesFEM,1)+1)-1;
                    tempIndexU(2*tempi)   = 2*(size(coordinatesFEM,1)+1);
                end
            end
            
            [ptOfxAll,ptOfyAll] = ndgrid(pt1x:pt3x,pt1y:pt3y);
            ptOfxAll = ptOfxAll(:); ptOfyAll = ptOfyAll(:); 
            
            %for ptOfx = pt1x:pt3x
            %    for ptOfy = pt1y:pt3y
            for tempjj = 1:length(ptOfxAll)
                
                ptOfx = ptOfxAll(tempjj); ptOfy = ptOfyAll(tempjj);
                    
                    % ------ Calculate ksi and eta ------
                    ksi = ksiMat(tempjj); eta = etaMat(tempjj);
                    %ksi = l(1)*ptOfx*ptOfy + l(2)*ptOfx + l(3)*ptOfy + l(4) ;
                    %eta = m(1)*ptOfx*ptOfy + m(2)*ptOfx + m(3)*ptOfy + m(4) ;
                    
                    % ------ Calculate N ------
                    deltaPt5 = 0; deltaPt6 = 0; deltaPt7 = 0; deltaPt8 = 0;
                    if elementsFEM(j,5) ~= 0, deltaPt5 = 1; end
                    if elementsFEM(j,6) ~= 0, deltaPt6 = 1; end
                    if elementsFEM(j,7) ~= 0, deltaPt7 = 1; end
                    if elementsFEM(j,8) ~= 0, deltaPt8 = 1; end
                    
                    N5 = deltaPt5*0.5*(1+ksi)*(1-abs(eta));
                    N6 = deltaPt6*0.5*(1+eta)*(1-abs(ksi));
                    N7 = deltaPt7*0.5*(1-ksi)*(1-abs(eta));
                    N8 = deltaPt8*0.5*(1-eta)*(1-abs(ksi));
                    
                    N1 = (1-ksi)*(1-eta)*0.25 - 0.5*(N7+N8);
                    N2 = (1+ksi)*(1-eta)*0.25 - 0.5*(N8+N5);
                    N3 = (1+ksi)*(1+eta)*0.25 - 0.5*(N5+N6);
                    N4 = (1-ksi)*(1+eta)*0.25 - 0.5*(N6+N7);
                    
                    N = [N1 0 N2 0 N3 0 N4 0 N5 0 N6 0 N7 0 N8 0;
                         0 N1 0 N2 0 N3 0 N4 0 N5 0 N6 0 N7 0 N8];
                    
                    % J11 = funDN1Dksi(ksi,eta,deltaPt7,deltaPt8)*pt1x + funDN2Dksi(ksi,eta,deltaPt8,deltaPt5)*pt2x + ...
                    %       funDN3Dksi(ksi,eta,deltaPt5,deltaPt6)*pt3x + funDN4Dksi(ksi,eta,deltaPt6,deltaPt7)*pt4x + ...
                    %       funDN5Dksi(ksi,eta,deltaPt5)*pt5x + funDN6Dksi(ksi,eta,deltaPt6)*pt6x + ...
                    %       funDN7Dksi(ksi,eta,deltaPt7)*pt7x + funDN8Dksi(ksi,eta,deltaPt8)*pt8x;
                    % J12 = funDN1Dksi(ksi,eta,deltaPt7,deltaPt8)*pt1y + funDN2Dksi(ksi,eta,deltaPt8,deltaPt5)*pt2y + ...
                    %       funDN3Dksi(ksi,eta,deltaPt5,deltaPt6)*pt3y + funDN4Dksi(ksi,eta,deltaPt6,deltaPt7)*pt4y + ...
                    %       funDN5Dksi(ksi,eta,deltaPt5)*pt5y + funDN6Dksi(ksi,eta,deltaPt6)*pt6y + ...
                    %       funDN7Dksi(ksi,eta,deltaPt7)*pt7y + funDN8Dksi(ksi,eta,deltaPt8)*pt8y;
                    % J21 = funDN1Deta(ksi,eta,deltaPt7,deltaPt8)*pt1x + funDN2Deta(ksi,eta,deltaPt8,deltaPt5)*pt2x + ...
                    %       funDN3Deta(ksi,eta,deltaPt5,deltaPt6)*pt3x + funDN4Deta(ksi,eta,deltaPt6,deltaPt7)*pt4x + ...
                    %       funDN5Deta(ksi,eta,deltaPt5)*pt5x + funDN6Deta(ksi,eta,deltaPt6)*pt6x + ...
                    %       funDN7Deta(ksi,eta,deltaPt7)*pt7x + funDN8Deta(ksi,eta,deltaPt8)*pt8x;
                    % J22 = funDN1Deta(ksi,eta,deltaPt7,deltaPt8)*pt1y + funDN2Deta(ksi,eta,deltaPt8,deltaPt5)*pt2y + ...
                    %       funDN3Deta(ksi,eta,deltaPt5,deltaPt6)*pt3y + funDN4Deta(ksi,eta,deltaPt6,deltaPt7)*pt4y + ...
                    %       funDN5Deta(ksi,eta,deltaPt5)*pt5y + funDN6Deta(ksi,eta,deltaPt6)*pt6y + ...
                    %       funDN7Deta(ksi,eta,deltaPt7)*pt7y + funDN8Deta(ksi,eta,deltaPt8)*pt8y;
                    % 
                    % J = [J11 J12; J21 J22];
                    
                    J = [funDN1Dksi(ksi,eta,deltaPt7,deltaPt8) funDN2Dksi(ksi,eta,deltaPt8,deltaPt5) ...
                          funDN3Dksi(ksi,eta,deltaPt5,deltaPt6) funDN4Dksi(ksi,eta,deltaPt6,deltaPt7) ...
                          funDN5Dksi(ksi,eta,deltaPt5) funDN6Dksi(ksi,eta,deltaPt6) ...
                          funDN7Dksi(ksi,eta,deltaPt7) funDN8Dksi(ksi,eta,deltaPt8);
                          funDN1Deta(ksi,eta,deltaPt7,deltaPt8) funDN2Deta(ksi,eta,deltaPt8,deltaPt5) ...
                          funDN3Deta(ksi,eta,deltaPt5,deltaPt6) funDN4Deta(ksi,eta,deltaPt6,deltaPt7) ...
                          funDN5Deta(ksi,eta,deltaPt5) funDN6Deta(ksi,eta,deltaPt6) ...
                          funDN7Deta(ksi,eta,deltaPt7) funDN8Deta(ksi,eta,deltaPt8)] * ...
                        [pt1x,pt1y;pt2x,pt2y;pt3x,pt3y;pt4x,pt4y;pt5x,pt5y;pt6x,pt6y;pt7x,pt7y;pt8x,pt8y];
                     
                    Jacobian = det(J);
                    InvJ = 1/Jacobian*[J(2,2) -J(1,2); -J(2,1) J(1,1)];
                    
                    DN = [InvJ zeros(2,2); zeros(2,2) InvJ] * ...
                        [funDN1Dksi(ksi,eta,deltaPt7,deltaPt8) 0 funDN2Dksi(ksi,eta,deltaPt8,deltaPt5) 0 ...
                         funDN3Dksi(ksi,eta,deltaPt5,deltaPt6) 0 funDN4Dksi(ksi,eta,deltaPt6,deltaPt7) 0 ...
                         funDN5Dksi(ksi,eta,deltaPt5) 0 funDN6Dksi(ksi,eta,deltaPt6) 0 ...
                         funDN7Dksi(ksi,eta,deltaPt7) 0 funDN8Dksi(ksi,eta,deltaPt8) 0;
                         funDN1Deta(ksi,eta,deltaPt7,deltaPt8) 0 funDN2Deta(ksi,eta,deltaPt8,deltaPt5) 0 ...
                         funDN3Deta(ksi,eta,deltaPt5,deltaPt6) 0 funDN4Deta(ksi,eta,deltaPt6,deltaPt7) 0 ...
                         funDN5Deta(ksi,eta,deltaPt5) 0 funDN6Deta(ksi,eta,deltaPt6) 0 ...
                         funDN7Deta(ksi,eta,deltaPt7) 0 funDN8Deta(ksi,eta,deltaPt8) 0;
                         0 funDN1Dksi(ksi,eta,deltaPt7,deltaPt8) 0 funDN2Dksi(ksi,eta,deltaPt8,deltaPt5) ...
                         0 funDN3Dksi(ksi,eta,deltaPt5,deltaPt6) 0 funDN4Dksi(ksi,eta,deltaPt6,deltaPt7) ...
                         0 funDN5Dksi(ksi,eta,deltaPt5) 0 funDN6Dksi(ksi,eta,deltaPt6) ...
                         0 funDN7Dksi(ksi,eta,deltaPt7) 0 funDN8Dksi(ksi,eta,deltaPt8);
                         0 funDN1Deta(ksi,eta,deltaPt7,deltaPt8) 0 funDN2Deta(ksi,eta,deltaPt8,deltaPt5) ...
                         0 funDN3Deta(ksi,eta,deltaPt5,deltaPt6) 0 funDN4Deta(ksi,eta,deltaPt6,deltaPt7) ...
                         0 funDN5Deta(ksi,eta,deltaPt5) 0 funDN6Deta(ksi,eta,deltaPt6) ...
                         0 funDN7Deta(ksi,eta,deltaPt7) 0 funDN8Deta(ksi,eta,deltaPt8);];
                    
                    % ------ Here approximate Dg(x+u)=Df(x) ------
                    DfEle = [DfDx(ptOfx-DfDxStartx, ptOfy-DfDxStarty); 
                             DfDy(ptOfx-DfDxStartx, ptOfy-DfDxStarty)];
                    
                      % ------ Only assemble stiffness in the first step ------
                    if (stepwithinwhile==1)
                        %A(temp,temp) =  A(temp,temp) + (N'*Df)*(N'*Df)' + alpha*(DN')*DN ;
                        tempA = tempA + (N'*DfEle)*(N'*DfEle)' + alpha*(DN')*DN;
                    end
                         
                    temp1 = [ptOfx;ptOfy] + N*U(tempIndexU);
                    temp2 = ((f(ptOfx ,ptOfy) - fungInterpolation_g(temp1(1), temp1(2), g(floor(temp1(1))-1:floor(temp1(1))+2, floor(temp1(2))-1:floor(temp1(2))+2))) * (N'*DfEle));
                    
                    % Localb(j,:) = reshape(tempb,1,2*size(coordinatesFEM,1));
                    % OLD VERSION: tempb = ((f(ptOfx ,ptOfy) - fungInterpolation_g(temp1(1), temp1(2), g)) * (N'*Df));
                    
                    % b(temp) = b(temp) + tempb;
                    %b(temp) = b(temp) + tempb - (alpha*(DN')*DN)*U(temp);
                    tempb = tempb + temp2 - (alpha*(DN')*DN)*U(tempIndexU);
                     
                    % for brow = 1:2*size(coordinatesFEM,1)
                    %    nFEMbVector = nFEMbVector + 1;
                    %    bVectorI(nFEMbVector) = temp(brow);
                    %    bVectorK(nFEMbVector) = tempb(brow);
                    % end
                    
                %end
            end
            
            
            % =====================================================================
            % --- To store A_ele for each element ---
            if (stepwithinwhile==1)
                % [IndexAXX,IndexAYY] = ndgrid(tempIndexU,tempIndexU);
                % INDEXAI = [INDEXAI;IndexAXX(:)]; INDEXAJ = [INDEXAJ;IndexAYY(:)]; INDEXAVAL = [INDEXAVAL;tempA(:)]; %INDEXAREG = [INDEXAREG;tempAreg(:)];
                [IndexAXX,IndexAYY] = ndgrid(tempIndexU,tempIndexU);
                INDEXAIpar{j}=IndexAXX(:); INDEXAJpar{j}=IndexAYY(:); INDEXAVALpar{j}=tempA(:);
            end
            %INDEXBI = [INDEXBI;tempIndexU(:)]; INDEXBVAL = [INDEXBVAL;tempb(:)];
            INDEXBIpar{j}=tempIndexU(:); INDEXBVALpar{j}=tempb(:);   
            
            %    A_ele(j,1:256) = reshape(A(temp,temp),1,256);
            %    % ------ No need to modify A matrix ---------
            %    % if elementsFEM(j,5) == (size(coordinatesFEM,1)+1)
            %    %     A_ele(j,16*8+1:16*10) = zeros(1,32);
            %    % end
            %    % if elementsFEM(j,6) == (size(coordinatesFEM,1)+1)
            %    %     A_ele(j,16*10+1:16*12) = zeros(1,32);
            %    % end
            %    % if elementsFEM(j,7) == (size(coordinatesFEM,1)+1)
            %    %     A_ele(j,16*12+1:16*14) = zeros(1,32);
            %    % end
            %    % if elementsFEM(j,8) == (size(coordinatesFEM,1)+1)
            %    %     A_ele(j,16*14+1:16*16) = zeros(1,32);
            %    % end
            %    % -------------------------------------------
            %
            % end
            % 
            % --- To use Navier-Lame elasticity regularization instead ---
            % b(temp) = b(temp) + (mu*Matrix1Update + (mu+lamda)*Matrix2Update)*U(temp);
            % =====================================================================
            
        end
    end
    
    close(hbar);
    
     
    
    % ================= UNUSED =============
    % for j = 1: size(elementsFEM,1)
    %    temp = Localindex(j,:);
    %    A(temp,temp) =  A(temp,temp) + LocalA(j,:,:);
    %    b(temp) = b(temp) + reshape(Localb(j,:),2*size(coordinatesFEM,1),1);
    % end
    
    % ================= UNUSED =============
    % if (stepwithinwhile==1)
    %    A = sparse(StiffnessMatrixI,StiffnessMatrixJ,StiffnessMatrixK);
    % end
    % b = sparse(bVectorI, ones(size(bVectorI,1),1), bVectorK);
    
    
    % ================= UNUSED =============
    % [W,LinearSystemSolveFlag,~,~,~]  = pcg(A,b+2*lamda*LK'*Lb,1e-8,100,L,L');
    % [W,LinearSystemSolveFlag,~,~,~]  = pcg(A,b,1e-8,100,L,L');
    % if (LinearSystemSolveFlag~=0)
    %     disp('something wrong with pcg solver');
    % end
    
    % Largest elementsFEM index is located at the end; and "unique" func in order
%     if stepwithinwhile == 1
%         A = sparse(INDEXAI,INDEXAJ,INDEXAVAL,FEMSize,FEMSize);
%     end
%     b = sparse(INDEXBI,ones(length(INDEXBI),1),INDEXBVAL,FEMSize,1);
if stepwithinwhile == 1
    A = sparse(FEMSize, FEMSize);
    for eleInd = 1:size(elementsFEM,1)
        A = A + sparse(INDEXAIpar{eleInd}, INDEXAJpar{eleInd}, INDEXAVALpar{eleInd}, FEMSize,FEMSize) ;
        
    end
end
b = sparse(FEMSize,1);
for eleInd = 1:size(elementsFEM,1)
    
    b = b + sparse(INDEXBIpar{eleInd},ones(length(INDEXBIpar{eleInd}),1),INDEXBVALpar{eleInd}, FEMSize,1) ;
end

    
    %AMatrixRegularized = (A+1e-3*max(diag(A))*eye(2*size(coordinatesFEM,1)+2,2*size(coordinatesFEM,1)+2));
    
    coordsIndexInvolved = unique(elementsFEM);
    UIndexInvolved = [coordsIndexInvolved(2:end);coordsIndexInvolved(2:end)];
    % Not including the first 0-th entry
    for tempi = 1:(size(coordsIndexInvolved,1)-1)
        UIndexInvolved(2*tempi-1:2*tempi) = [2*coordsIndexInvolved(tempi+1)-1; 2*coordsIndexInvolved(tempi+1)];
    end
     
    W = sparse(2*size(coordinatesFEM,1),1);
    W(2*unique(dirichlet)) = 0;
    W(2*unique(dirichlet)-1) = 0;
    
    dirichlettemp = [2*dirichlet; 2*dirichlet-1];
    FreeNodes = setdiff(UIndexInvolved,unique(dirichlettemp));
    
    W(FreeNodes) = A(FreeNodes,FreeNodes)\b(FreeNodes);
      
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
    
    
     
    if size(UIndexInvolved,1) == length(W)
        % if (normW < 1e-6 || ((normW/normWOld > 0.9) && (normW/normWOld < 1)))
        %     U(UIndexInvolved(1:end)) = U(UIndexInvolved(1:end)) + W;
        %     break;   
        % elseif ((normW >= 1e-6) && (normW/normWOld < 1.5))
        %     U(UIndexInvolved(1:end)) = U(UIndexInvolved(1:end)) + W;
        % else
        %     warning('Get diverged in Global_ICGN!!!')
        %     break; 
        % end
        if (normW < tol) || ((normW/normWOld > 1-tol) && (normW/normWOld < 1))
            U(UIndexInvolved(1:end)) = U(UIndexInvolved(1:end)) + W;
            break;
        elseif  (normW >= tol && normW < 1/tol)  
            U(UIndexInvolved(1:end)) = U(UIndexInvolved(1:end)) + W;
        else
            warning('Get diverged in Global_ICGN!!!')
            break;
        end
    else
        warning('Index UIndexInvolved isnt matching with length of W vector!!!')
    end
    
end
    
 

U = U(1:end-2);

end


%% ========= subroutines for  FEM Q4 shape function derivatives ========

function DN5Dksi=funDN5Dksi(ksi,eta,deltaPt5)
DN5Dksi = deltaPt5*0.5*(1-abs(eta)) ;
end
function DN5Deta=funDN5Deta(ksi,eta,deltaPt5)
DN5Deta = deltaPt5*0.5*(1+ksi)*sign(-eta);
end
function DN6Dksi=funDN6Dksi(ksi,eta,deltaPt6)
DN6Dksi = deltaPt6*0.5*(1+eta)*sign(-ksi);
end
function DN6Deta=funDN6Deta(ksi,eta,deltaPt6)
DN6Deta = deltaPt6*0.5*(1-abs(ksi));
end
function DN7Dksi=funDN7Dksi(ksi,eta,deltaPt7)
DN7Dksi = deltaPt7*0.5*(-1)*(1-abs(eta));
end
function DN7Deta=funDN7Deta(ksi,eta,deltaPt7)
DN7Deta = deltaPt7*0.5*(1-ksi)*sign(-eta);
end
function DN8Dksi=funDN8Dksi(ksi,eta,deltaPt8)
DN8Dksi = deltaPt8*0.5*(1-eta)*sign(-ksi);
end
function DN8Deta=funDN8Deta(ksi,eta,deltaPt8)
DN8Deta = deltaPt8*0.5*(-1)*(1-abs(ksi));
end
function DN1Dksi = funDN1Dksi(ksi,eta,deltaPt7,deltaPt8)
DN1Dksi = -0.25*(1-eta)-0.5*((deltaPt7*0.5*(-1)*(1-abs(eta)))+deltaPt8*0.5*(1-eta)*sign(-ksi));
end
function DN1Deta = funDN1Deta(ksi,eta,deltaPt7,deltaPt8)
DN1Deta = -0.25*(1-ksi)-0.5*(deltaPt7*0.5*(1-ksi)*sign(-eta)+deltaPt8*0.5*(-1)*(1-abs(ksi)));
end
function DN2Dksi = funDN2Dksi(ksi,eta,deltaPt8,deltaPt5)
DN2Dksi = 0.25*(1-eta)-0.5*(deltaPt8*0.5*(1-eta)*sign(-ksi)+deltaPt5*0.5*(1-abs(eta)));
end
function DN2Deta = funDN2Deta(ksi,eta,deltaPt8,deltaPt5)
DN2Deta = -0.25*(1+ksi)-0.5*(deltaPt8*0.5*(-1)*(1-abs(ksi))+deltaPt5*0.5*(1+ksi)*sign(-eta));
end
function DN3Dksi = funDN3Dksi(ksi,eta,deltaPt5,deltaPt6)
DN3Dksi = 0.25*(1+eta)-0.5*(deltaPt5*0.5*(1-abs(eta))+deltaPt6*0.5*(1+eta)*sign(-ksi));
end
function DN3Deta = funDN3Deta(ksi,eta,deltaPt5,deltaPt6)
DN3Deta = 0.25*(1+ksi)-0.5*(deltaPt5*0.5*(1+ksi)*sign(-eta)+deltaPt6*0.5*(1-abs(ksi))); 
end
function DN4Dksi = funDN4Dksi(ksi,eta,deltaPt6,deltaPt7)
DN4Dksi = -0.25*(1+eta)-0.5*(deltaPt6*0.5*(1+eta)*sign(-ksi)+deltaPt7*0.5*(-1)*(1-abs(eta)));
end
function DN4Deta = funDN4Deta(ksi,eta,deltaPt6,deltaPt7)
DN4Deta = 0.25*(1-ksi)-0.5*(deltaPt6*0.5*(1-abs(ksi))+deltaPt7*0.5*(1-ksi)*sign(-eta));
end

