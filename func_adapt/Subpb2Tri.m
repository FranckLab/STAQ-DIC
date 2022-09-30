%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function AL-DIC Subproblem 2                             %
% Object: to find deformation field using global methods   %
% Author: Jin Yang                                         %
% Last date modified: 2018.03                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Uhat,EnrHAndTipEleIndex,EnrTipEleIndex] = Subpb2Tri(coordinatesFEM, elementsFEM, dirichlet, neumann, beta, mu, ...
    U, F, W, v, alpha, clusterNo)

% ====== Set boundary conditions ======
% % No Dirichlet boundary conditions
% dirichlet = []; % Set boundary condition by default that dirichlet = [];
% Set Neumann boundary conditions using input "F" gradient tensors.
winstepsize = abs(coordinatesFEM(1,1) - coordinatesFEM(2,1));

% ====== Info of CrackTip path-line ======
% k_tip = -0.5*(CrackPath1(1)/CrackPath1(2)+CrackPath2(1)/CrackPath2(2));
% CrackPathCen = [k_tip/(CrackTip(2)-k_tip*CrackTip(1)), -1/(CrackTip(2)-k_tip*CrackTip(1))];

try
    % ------ Pre-Subpb2 deal with cracks ------
    % EnrHAndTipEleIndex is the collection of index for Heaviside function basis;
    % EnrTipEleIndex is the collection of index for Crack tip singular basis;
    if CrackOrNot == 1
        [EnrHAndTipEleIndex,EnrTipEleIndex] = EnrEleWithFrac(coordinatesFEM,elementsFEM,CrackPath1,CrackPath2,CrackTip);
    else
        EnrHAndTipEleIndex = []; EnrTipEleIndex = [];
    end
    
    if CrackOrNot == 0
        FEMSize = size(coordinatesFEM,1);
        U = U; F = F; W = 0*W; v = 0*v;
        EnrHAndTipEleIndex=[]; EnrTipEleIndex=[];
    elseif CrackTipOrNot == 1
        FEMSize = size(coordinatesFEM,1) + size(EnrHAndTipEleIndex,1) + 4*size(EnrTipEleIndex,1);
        W = 0*F; v = 0*U;
        F = [F;zeros(4*size(EnrHAndTipEleIndex,1),1);zeros(4*4*size(EnrTipEleIndex,1),1)];
        U = [U;zeros(2*size(EnrHAndTipEleIndex,1),1);zeros(2*4*size(EnrTipEleIndex,1),1)];
        W = [W;zeros(4*size(EnrHAndTipEleIndex,1),1);zeros(4*4*size(EnrTipEleIndex,1),1)];
        v = [v;zeros(2*size(EnrHAndTipEleIndex,1),1);zeros(2*4*size(EnrTipEleIndex,1),1)];
    else
        FEMSize = size(coordinatesFEM,1) + size(EnrHAndTipEleIndex,1);
        W = 0*F; v = 0*U;
        F = [F;zeros(4*size(EnrHAndTipEleIndex,1),1) ];
        U = [U;zeros(2*size(EnrHAndTipEleIndex,1),1) ];
        W = [W;zeros(4*size(EnrHAndTipEleIndex,1),1) ];
        v = [v;zeros(2*size(EnrHAndTipEleIndex,1),1) ];
    end
catch
    FEMSize = size(coordinatesFEM,1);
    W = 0*W; v = 0*v;
    EnrHAndTipEleIndex=[]; EnrTipEleIndex=[];
    CrackOrNot = 0; CrackTipOrNot = 0; % No crack, no crack tip
end
%%%%%%%%%%% !!!!!!!!!!!!! %%%%%%%%%%%%%%

% ====== Initialize variables ======
Uhat = U; UMinusv = U-v; FMinusW = F-W;

% ====== Initialize A matrix and b vector ======
A = sparse(2*FEMSize, 2*FEMSize); % 8 is because I put all the zeros in elementsFEM to the end, and in NC function there are 8 entries

% ====== Gaussian quadrature parameter ======
% ------ 4*4 Gaussian points ------
% gqpt1 = 0.339981; gqpt2 = -0.339981; gqpt3 = 0.861136; gqpt4 = -0.861136;
% gqwt1 = 0.652145; gqwt2 = 0.652145; gqwt3 = 0.347855; gqwt4 = 0.347855;
% gqpt = [gqpt1,gqpt2,gqpt3,gqpt4]; gqwt = [gqwt1,gqwt2,gqwt3,gqwt4];
% ------ 5*5 Gaussian points ------
% gqpt1 = 0; gqpt2 = 0.538469; gqpt3 = -0.538469; gqpt4 = 0.90618; gqpt5 = -0.90618;
% gqwt1 = 0.568889; gqwt2 = 0.478629; gqwt3 = 0.478629; gqwt4 = 0.236927; gqwt5 = 0.236927;
% gqpt = [gqpt1,gqpt2,gqpt3,gqpt4,gqpt5]; gqwt = [gqwt1,gqwt2,gqwt3,gqwt4,gqwt5];
% ------ 3*3 Gaussian points ------
% gqpt1 = 0; gqpt2 = 0.57735; gqpt3 = -0.57735; gqpt = [gqpt1,gqpt2,gqpt3];
% gqwt1 = 1; gqwt2 = 1; gqwt3 = 1; gqwt = [gqwt1,gqwt2,gqwt3];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ====== Initialize Global FEM solver ======
ConvergeOrNot = 0; IterStep = 0;

while ConvergeOrNot < 0.5 && IterStep < 1
    
    % ====== Initialize A matrix and b vector ======
    clear b; b = sparse(2*FEMSize,1); IterStep = IterStep+1;
    hbar = waitbar(0,['Global step']);
    % ============= Each element, assemble stiffness matrix ============
    for j = 1:size(elementsFEM,1)
        
        waitbar(j/size(elementsFEM,1));
        EleCrackTipOrNot = 0; EleCrackHOrNot = 0;
        
        % ------ Find three corner points ------
        point1x = coordinatesFEM(elementsFEM(j,1),1); point1y = coordinatesFEM(elementsFEM(j,1),2);
        point2x = coordinatesFEM(elementsFEM(j,2),1); point2y = coordinatesFEM(elementsFEM(j,2),2);
        point3x = coordinatesFEM(elementsFEM(j,3),1); point3y = coordinatesFEM(elementsFEM(j,3),2);
        
        % ------ Compute triangle area ------
        TriArea =  det([1 point1x point1y; 1 point2x point2y; 1 point3x point3y]);
        
        % ------ Compute DN matrix for CST element ------
        funDN1x = 1/(2*TriArea)*(point2y-point3y); funDN1y = 1/(2*TriArea)*(point3x-point2x);
        funDN2x = 1/(2*TriArea)*(point3y-point1y); funDN2y = 1/(2*TriArea)*(point1x-point3x);
        funDN3x = 1/(2*TriArea)*(point1y-point2y); funDN3y = 1/(2*TriArea)*(point2x-point1x);
        
        DN = [funDN1x,0,funDN2x,0,funDN3x,0; funDN1y,0,funDN2y,0,funDN3y,0;
            0 funDN1x,0,funDN2x,0,funDN3x ; 0 funDN1y,0,funDN2y,0,funDN3y ];
        
        % ------ Find the element nodal indices ------
        temp = [2*elementsFEM(j,1)-1, 2*elementsFEM(j,1), 2*elementsFEM(j,2)-1, 2*elementsFEM(j,2), 2*elementsFEM(j,3)-1, 2*elementsFEM(j,3)];
        
        tempF = [4*elementsFEM(j,1)-3 4*elementsFEM(j,1)-1 4*elementsFEM(j,1)-2 4*elementsFEM(j,1);
            4*elementsFEM(j,1)-3 4*elementsFEM(j,1)-1 4*elementsFEM(j,1)-2 4*elementsFEM(j,1);
            4*elementsFEM(j,2)-3 4*elementsFEM(j,2)-1 4*elementsFEM(j,2)-2 4*elementsFEM(j,2);
            4*elementsFEM(j,2)-3 4*elementsFEM(j,2)-1 4*elementsFEM(j,2)-2 4*elementsFEM(j,2);
            4*elementsFEM(j,3)-3 4*elementsFEM(j,3)-1 4*elementsFEM(j,3)-2 4*elementsFEM(j,3);
            4*elementsFEM(j,3)-3 4*elementsFEM(j,3)-1 4*elementsFEM(j,3)-2 4*elementsFEM(j,3)];
        
        % ------ Set Gauss points ------
        pointOfx = 1/3*(point1x+point2x+point3x);
        pointOfy = 1/3*(point1y+point2y+point3y);
        
        % Judge point is inside triangle or not
        % pointInTriangleOrNot = funPointInTriangleCheck(point1x,point1y,point2x,point2y,point3x,point3y,pointOfx,pointOfy);
        
        % if pointInTriangleOrNot == 1
        
        % ------ Calculate N ------
        N1 = det([1 pointOfx pointOfy; 1 point2x point2y; 1 point3x point3y])/TriArea;
        N2 = det([1 pointOfx pointOfy; 1 point3x point3y; 1 point1x point1y])/TriArea;
        N3 = det([1 pointOfx pointOfy; 1 point1x point1y; 1 point2x point2y])/TriArea;
        N = [N1, 0, N2, 0, N3, 0; 0, N1, 0, N2, 0, N3];
        NDiag = diag([N1,N1,N2,N2,N3,N3]);
        
        % ------ Construct A matrix ------
        %    if IterStep == 1
        A(temp,temp) = A(temp,temp) + TriArea*( (beta+alpha)*(DN')*(DN)+ mu*(N')*N ) ;
        %    end
        
        % ------ Construct b vector ------
        b(temp) = b(temp) + TriArea*( -beta*diag((DN')*(FMinusW(tempF)')) + mu*NDiag*(UMinusv(temp)) +(alpha)*(DN')*DN*U(temp) );
        
        % end
        
        
        
        
    end
    close(hbar);
    
    %  if IterStep == 1
    AMatrixRegularized = A + 1e-3*max(diag(A))*speye(size(A,1),size(A,2));
    % AMatrixRegularized = A; tempVal = 1e-3*max(diag(A));
    % for tempi = 1:size(A,1)
    %     AMatrixRegularized(tempi,tempi) = AMatrixRegularized(tempi,tempi) + tempVal;
    % send
    %  end
    
    %     coordsIndexInvolved = unique(elementsFEM); % Need modification for triangle elementsFEM
    %     if CrackOrNot == 1
    %         coordsIndexInvolvedH  = EnrHAndTipEleIndex(:,2);
    %         if CrackTipOrNot == 0
    %             coordsIndexInvolvedC1 = EnrTipEleIndex(:,2); coordsIndexInvolvedC2 = EnrTipEleIndex(:,2)+1;
    %             coordsIndexInvolvedC3 = EnrTipEleIndex(:,2)+2; coordsIndexInvolvedC4 = EnrTipEleIndex(:,2)+3;
    %             coordsIndexInvolved   = [coordsIndexInvolved;coordsIndexInvolvedH;coordsIndexInvolvedC1;coordsIndexInvolvedC2;coordsIndexInvolvedC3;coordsIndexInvolvedC4];
    %         else
    %             coordsIndexInvolved   = [coordsIndexInvolved;coordsIndexInvolvedH];
    %         end
    %     end
    %     UIndexInvolved = zeros(2*(length(coordsIndexInvolved)),1);
    %     % Not including the first 0-th entry
    %     for tempi = 1:(size(coordsIndexInvolved,1))
    %         UIndexInvolved(2*tempi-1:2*tempi) = [2*coordsIndexInvolved(tempi)-1; 2*coordsIndexInvolved(tempi)];
    %     end
    dirichlettemp = [2*dirichlet; 2*dirichlet-1];
    FreeNodes = setdiff(1:2*size(coordinatesFEM,1),unique(dirichlettemp));
    %     % ========= Set Dirichlet and Neumann boundary conditions =========
    %     if isempty(dirichlet) ~= 1
    %         dirichlettemp = [2*dirichlet(:); 2*dirichlet(:)-1];
    %     else
    %         dirichlettemp = [];
    %     end
    %     if isempty(neumann) ~= 1
    %         neumanntemp = [2*neumann(:,1); 2*neumann(:,1)-1; 2*neumann(:,2); 2*neumann(:,2)-1];
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
    UhatOld = Uhat;
    Uhat = sparse(2*FEMSize ,1);
    Uhat(2*unique(dirichlet)) = U(2*unique(dirichlet));
    Uhat(2*unique(dirichlet)-1) = U(2*unique(dirichlet)-1);
    b = b - AMatrixRegularized * Uhat;
    
    % [tempgridx,tempgridy] = meshgrid(1:FEMSize,1:FEMSize);
    % figure; mesh(tempgridx,tempgridy,AMatrixRegularized(1:FEMSize,1:FEMSize));
    
    % ========= Neumann boundary conditions ==========
    % Get boundary traction force using displacement input "U"
    BCForce = -1/(winstepsize)*F; % F is the input of Subpb2.m, which is the Subproblem 1 solved gradient def tensors;
    
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
    
    
    % ========= Solve FEM problem ===========
    Uhat(FreeNodes) = AMatrixRegularized(FreeNodes,FreeNodes) \ b(FreeNodes);
    UhatNew = Uhat(1:end );
    
    if norm(UhatNew-UhatOld)/sqrt(length(UhatOld)) < 1e-3
        ConvergeOrNot = 1;
    end
    
    Uhat = UhatNew;
    
end



%% ========= subroutines for  FEM shape function derivatives ========



