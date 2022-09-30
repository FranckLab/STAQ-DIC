function [Uhat] = Subpb2(DICmesh,beta,mu,U,F,udual,vdual,alpha,GaussPtOrder)
%FUNCTION [Uhat] = Subpb2(DICmesh,beta,mu,U,F,udual,vdual,alpha,GaussPtOrder)
% AL-DIC Subproblem 2 is solved over a uniform FE-mesh to find a globally
% kinematically compatible deformation field by finite element method.
% ----------------------------------------------
% 
%   INPUT: DICmesh             DIC FE Q4 mesh: coordinatesFEM, elementsFEM
%          GaussPtOrder        Gauss point order used in FE Q4 element, GaussPtOrder = 3 (default)
%          beta, mu            Two constant coefficients
%          U                   Disp vector: U = [Ux_node1, Uy_node1, ... , Ux_nodeN, Uy_nodeN]';
%          F                   Deformation gradient: F = [F11_node1, F21_node1, F12_node1, F22_node1, ... 
%                                                         ... , F11_nodeN, F21_nodeN, F12_nodeN, F22_nodeN]';
%          udual,vdual         Dual variables
%          alpha               Smoothness coefficient. Not needed here, i.e., alpha=0
%           
%   OUTPUT: Uhat               Solved globally kinematically compatible displacement field
%
% ----------------------------------------------
% Author: Jin Yang.  
% Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
% Last time updated: 2018.03, 2020.12 
% ==============================================  


%% Initialization
coordinatesFEM = DICmesh.coordinatesFEM;
elementsFEM = DICmesh.elementsFEM;
M = DICmesh.M;   N = DICmesh.N;
winstepsize = abs(coordinatesFEM(1,1) - coordinatesFEM(2,1));

DIM = 2;
NodesPerEle = 4;
FEMSize = size(coordinatesFEM,1);

% ====== Set boundary conditions ======
dirichlet = []; % Set boundary condition by default that dirichlet = [];
neumann = DICmesh.neumann; % Set Neumann boundary conditions using input "F" gradient tensors.
 

%% ====== Initialize FE-system ======
% Please ignore these codes. 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A = sparse(DIM*FEMSize,DIM*FEMSize);
% b = sparse(DIM*FEMSize,1);

% ====== Compute Div(F-W) ======
% FMinusW = F - W;
% FMinusW11 = reshape(FMinusW(1:4:end),M,N); FMinusW21 = reshape(FMinusW(2:4:end),M,N);
% FMinusW12 = reshape(FMinusW(3:4:end),M,N); FMinusW22 = reshape(FMinusW(4:4:end),M,N);

% imgradientMatrix = [-1/2 0 1/2]'; % Central finite difference operator
% DFMinusW11Dxtemp = imfilter(FMinusW11,imgradientMatrix);
% DFMinusW12Dytemp = imfilter(FMinusW12,imgradientMatrix');
% DFMinusW11Dx     = DFMinusW11Dxtemp(1:end, 1:end);
% DFMinusW12Dy     = DFMinusW12Dytemp(1:end, 1:end);
%
% DFMinusW21Dxtemp = imfilter(FMinusW21,imgradientMatrix);
% DFMinusW22Dytemp = imfilter(FMinusW22,imgradientMatrix');
% DFMinusW21Dx     = DFMinusW21Dxtemp(1:end, 1:end);
% DFMinusW22Dy     = DFMinusW22Dytemp(1:end, 1:end);
%
% DivFMinusW1temp = DFMinusW11Dx + DFMinusW12Dy;
% DivFMinusW2temp = DFMinusW21Dx + DFMinusW22Dy;
%
% DivFMinusWGlobal  = zeros(M*N*2,1);
% DivFMinusWGlobal(1:2:end) = reshape(DivFMinusW1temp,M*N,1);
% DivFMinusWGlobal(2:2:end) = reshape(DivFMinusW2temp,M*N,1);
%
% DFMinusW11DxStartx = 1; DFMinusW11DyStarty = 1;


%%
INDEXAI = []; INDEXAJ = []; INDEXAVAL = []; INDEXBI = []; INDEXBVAL = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ============= Each element, assemble stiffness matrix ============
for eleInd = 1:size(elementsFEM,1)
    
    % ------ Find four corner points ------
    pt1x = coordinatesFEM(elementsFEM(eleInd,1),1); pt1y = coordinatesFEM(elementsFEM(eleInd,1),2);
    pt2x = coordinatesFEM(elementsFEM(eleInd,2),1); pt2y = coordinatesFEM(elementsFEM(eleInd,2),2);
    pt3x = coordinatesFEM(elementsFEM(eleInd,3),1); pt3y = coordinatesFEM(elementsFEM(eleInd,3),2);
    pt4x = coordinatesFEM(elementsFEM(eleInd,4),1); pt4y = coordinatesFEM(elementsFEM(eleInd,4),2);
    
    % ------ Find the element nodal indices ------
    tempIndexU = 2*elementsFEM(eleInd,[1,1,2,2,3,3,4,4]);
    tempIndexU(1:2:end) = tempIndexU(1:2:end)-1;
    % ------ Find deformation gradient tensors indices ------
    tempIndexF = [4*elementsFEM(eleInd,[1,1,2,2,3,3,4,4])'-3, ...
                  4*elementsFEM(eleInd,[1,1,2,2,3,3,4,4])'-1, ...
                  4*elementsFEM(eleInd,[1,1,2,2,3,3,4,4])'-2, ...
                  4*elementsFEM(eleInd,[1,1,2,2,3,3,4,4])'];
                  
    % ------ Compute DivFMinusW ------
    FMinusW = F(tempIndexF) - udual(tempIndexF);
     
    % ------ Compute UMinusV ------
    UMinusV = U(tempIndexU) - vdual(tempIndexU);
     
    % ------ Nine Gauss integral points ------
    ksietaList = [-sqrt(3/5), 0, sqrt(3/5)];
    ksietaWeightList = [5/9, 8/9, 5/9];
    
    tempA = zeros(DIM*NodesPerEle,DIM*NodesPerEle);  tempb = zeros(DIM*NodesPerEle,1);
    
    for tempi = 1:3
        for tempj = 1:3
            ksi = ksietaList(tempi);
            eta = ksietaList(tempj);
            weightksi = ksietaWeightList(tempi);
            weighteta = ksietaWeightList(tempj);
            
            % ------ Calculate N matrix ------
            N1 = (1-ksi)*(1-eta)*0.25;
            N2 = (1+ksi)*(1-eta)*0.25;
            N3 = (1+ksi)*(1+eta)*0.25;
            N4 = (1-ksi)*(1+eta)*0.25;
            N = [N1 0 N2 0 N3 0 N4 0;
                0 N1 0 N2 0 N3 0 N4];
            NDiag = diag([N1,N1,N2,N2,N3,N3,N4,N4]);
            
            % ------ Calculate Jacobian matrix ------
            J11 = funDN1Dksi(ksi,eta)*pt1x + funDN2Dksi(ksi,eta)*pt2x + ...
                  funDN3Dksi(ksi,eta)*pt3x + funDN4Dksi(ksi,eta)*pt4x;
            J12 = funDN1Dksi(ksi,eta)*pt1y + funDN2Dksi(ksi,eta)*pt2y + ...
                  funDN3Dksi(ksi,eta)*pt3y + funDN4Dksi(ksi,eta)*pt4y;
            J21 = funDN1Deta(ksi,eta)*pt1x + funDN2Deta(ksi,eta)*pt2x + ...
                  funDN3Deta(ksi,eta)*pt3x + funDN4Deta(ksi,eta)*pt4x;
            J22 = funDN1Deta(ksi,eta)*pt1y + funDN2Deta(ksi,eta)*pt2y + ...
                  funDN3Deta(ksi,eta)*pt3y + funDN4Deta(ksi,eta)*pt4y;
            J = [J11 J12; J21 J22];
            
            Jacobian = det(J);
            InvJ = 1/Jacobian*[J22 -J12; -J21 J11];
            
            % ------ Compute DN matrix ------
            DN = [InvJ zeros(2,2); zeros(2,2) InvJ] * ...
                [funDN1Dksi(ksi,eta) 0 funDN2Dksi(ksi,eta) 0 funDN3Dksi(ksi,eta) 0 funDN4Dksi(ksi,eta) 0;
                funDN1Deta(ksi,eta) 0 funDN2Deta(ksi,eta) 0 funDN3Deta(ksi,eta) 0 funDN4Deta(ksi,eta) 0;
                0 funDN1Dksi(ksi,eta) 0 funDN2Dksi(ksi,eta) 0 funDN3Dksi(ksi,eta) 0 funDN4Dksi(ksi,eta);
                0 funDN1Deta(ksi,eta) 0 funDN2Deta(ksi,eta) 0 funDN3Deta(ksi,eta) 0 funDN4Deta(ksi,eta)];
            
            % ------ Construct A matrix ------
            %A(tempIndexU,tempIndexU) = A(tempIndexU,tempIndexU) + Jacobian*weightksi*weighteta * ( beta*(DN')*(DN) + mu*(N')*N + alpha*mu*(DN')*(DN) );
            tempA = tempA + Jacobian*weightksi*weighteta * ( beta*(DN')*(DN) + mu*(N')*N + alpha*mu*(DN')*(DN) );
            
            % ------ Construct b vector ------
            % b(temp) = b(temp) + Jacobian*weightksi*weighteta* ( -beta*NDiag*(DivFMinusW) + mu*NDiag*(UMinusV) );
            % b(tempIndexU) = b(tempIndexU) + Jacobian*weightksi*weighteta* ( beta*diag((DN')*(FMinusW')) + mu*(N')*N*(UMinusV) );
            tempb = tempb + Jacobian*weightksi*weighteta* ( beta*diag((DN')*(FMinusW')) + mu*(N')*N*(UMinusV) );
            
        end
    end
    
    % --- To tempA and tempb ---
    [IndexAXX,IndexAYY] = ndgrid(tempIndexU,tempIndexU);
    INDEXAI = [INDEXAI;IndexAXX(:)]; INDEXAJ = [INDEXAJ;IndexAYY(:)]; INDEXAVAL = [INDEXAVAL;tempA(:)];
    INDEXBI = [INDEXBI;tempIndexU(:)]; INDEXBVAL = [INDEXBVAL;tempb(:)];
    
end

% ========= Assemble A matrix and b vector =========
A = sparse(INDEXAI,INDEXAJ,INDEXAVAL, DIM*FEMSize, DIM*FEMSize);
b = sparse(INDEXBI, ones(length(INDEXBI),1),INDEXBVAL, DIM*FEMSize,1 );

% ========= Dirichlet boundary conditions ==========
Uhat = sparse(2*size(coordinatesFEM,1),1);
Uhat(2*unique(dirichlet)) = U(2*unique(dirichlet));
Uhat(2*unique(dirichlet)-1) = U(2*unique(dirichlet)-1);
b = b - A * Uhat;

dirichlettemp = [2*dirichlet; 2*dirichlet-1];
FreeNodes = setdiff(1:2*size(coordinatesFEM,1),unique(dirichlettemp));

% ========= Neumann boundary conditions ==========
% Get boundary traction force using displacement input "U"
BCForce = -1/(winstepsize)*F; % F is the input of Subpb2.m, which is the Subproblem 1 solved gradient def tensors;

% If F solved by Subproblem 1 is too noisy at boundary, try the following code:
% [BCForce] = -1/(winstepsize*2)*funGlobal_NodalStrainAvg(coordinatesFEM,elementsFEM,U,GaussPtOrder);
for tempj = 1:size(neumann,1)
    % f1 = F11*n1 + F12*n2
    b(2*neumann(tempj,1:2)-1) = b(2*neumann(tempj,1:2)-1) + 0.5*norm(coordinatesFEM(neumann(tempj,1),:)-coordinatesFEM(neumann(tempj,2),:)) ...
        *( ( BCForce(4*neumann(tempj,1:2)-3) * neumann(tempj,3) + BCForce(4*neumann(tempj,1:2)-1) * neumann(tempj,4) ) );
    % f2 = F21*n1 + F22*n2
    b(2*neumann(tempj,1:2))   = b(2*neumann(tempj,1:2)) + 0.5*norm(coordinatesFEM(neumann(tempj,1),:)-coordinatesFEM(neumann(tempj,2),:)) ...
        *( ( BCForce(4*neumann(tempj,1:2)-2) * neumann(tempj,3) + BCForce(4*neumann(tempj,1:2)) * neumann(tempj,4) ) ) ;
    
    %b(2*neumann(tempj,1:2)-1) = b(2*neumann(tempj,1:2)-1) + 0.5*1 ...
    %    *( ( BCForce(4*neumann(tempj,1:2)-3) * neumann(tempj,3) + BCForce(4*neumann(tempj,1:2)-1) * neumann(tempj,4) ) );
    %b(2*neumann(tempj,1:2))   = b(2*neumann(tempj,1:2)) + 0.5*1 ...
    %    *( ( BCForce(4*neumann(tempj,1:2)-2) * neumann(tempj,3) + BCForce(4*neumann(tempj,1:2)) * neumann(tempj,4) ) ) ;
    
end

Uhat(FreeNodes) = A(FreeNodes,FreeNodes) \ b(FreeNodes);

end


%% ========= subroutines for FEM Q4 shape function derivatives ========
function DN1Dksi=funDN1Dksi(ksi,eta)
    DN1Dksi = -(1-eta)/4 ;
end
function DN1Deta=funDN1Deta(ksi,eta)
    DN1Deta =  -(1-ksi)/4 ;
end
function DN2Dksi=funDN2Dksi(ksi,eta)
    DN2Dksi =  (1-eta)/4 ;
end
function DN2Deta=funDN2Deta(ksi,eta)
    DN2Deta =  -(1+ksi)/4 ;
end
function DN3Dksi=funDN3Dksi(ksi,eta)
    DN3Dksi = (1+eta)/4 ;
end
function DN3Deta=funDN3Deta(ksi,eta)
    DN3Deta =  (1+ksi)/4 ;
end
function DN4Dksi=funDN4Dksi(ksi,eta)
    DN4Dksi = -(1+eta)/4 ;
end
function DN4Deta=funDN4Deta(ksi,eta)
    DN4Deta = (1-ksi)/4 ;
end

