%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function AL-DIC Subproblem 2                             %
% Object: to find deformation field using global methods   %
% Author: Jin Yang                                         %
% Last date modified: 2018.03                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Uhat,EnrHAndTipEleIndex, EnrTipEleIndex] = Subpb2_adapt_para(coordinatesFEM, elementsFEM, dirichlet, neumann, beta, mu, ...
    USubpb1, FSubpb1, udual, vdual, alpha, alpha2, winstepsize, ...
    CrackOrNot,CrackPath1,CrackPath2,CrackTip,CrackTipOrNot,ALSolveStep,clusterNo)
 
% disp(['--- Start step ',num2str(ALSolveStep),' subproblem 2 ---']);
DIM = 2; NodesNumPerEle = 4; GaussPtOrder = 2; 

% ====== Info of CrackTip path-line ======
k_tip = -0.5*(CrackPath1(1)/CrackPath1(2)+CrackPath2(1)/CrackPath2(2));
CrackPathCen = [k_tip/(CrackTip(2)-k_tip*CrackTip(1)), -1/(CrackTip(2)-k_tip*CrackTip(1))];

% ------ Pre-Subpb2 deal with cracks ------
% EnrHAndTipEleIndex is the collection of index for Heaviside function basis;
% EnrTipEleIndex is the collection of index for Crack tip singular basis; 
% ------ ^^ I'm a ---- >< ---- 666 ---- cutting line ---------
if CrackOrNot == 1 
    [EnrHAndTipEleIndex,EnrTipEleIndex] = EnrEleWithFrac(coordinatesFEM,elementsFEM,CrackPath1,CrackPath2,CrackTip);
else
    EnrHAndTipEleIndex = []; EnrTipEleIndex = [];
end
if CrackOrNot == 0
    FEMSize = DIM*size(coordinatesFEM,1);
    U = USubpb1; F = FSubpb1; W = 0*udual; v = 0*vdual;
    EnrHAndTipEleIndex=[]; EnrTipEleIndex=[]; 
elseif CrackTipOrNot == 1
    FEMSize = DIM*(size(coordinatesFEM,1) + size(EnrHAndTipEleIndex,1) + 4*size(EnrTipEleIndex,1));
    udual = 0*FSubpb1; vdual = 0*USubpb1;
    F = [FSubpb1;zeros(4*size(EnrHAndTipEleIndex,1),1);zeros(4*4*size(EnrTipEleIndex,1),1)];
    U = [USubpb1;zeros(2*size(EnrHAndTipEleIndex,1),1);zeros(2*4*size(EnrTipEleIndex,1),1)];
    W = [udual;zeros(4*size(EnrHAndTipEleIndex,1),1);zeros(4*4*size(EnrTipEleIndex,1),1)];
    v = [vdual;zeros(2*size(EnrHAndTipEleIndex,1),1);zeros(2*4*size(EnrTipEleIndex,1),1)];
else
    FEMSize = DIM*(size(coordinatesFEM,1) + size(EnrHAndTipEleIndex,1));
    udual = 0*FSubpb1; vdual = 0*USubpb1;
    F = [FSubpb1;zeros(4*size(EnrHAndTipEleIndex,1),1) ];
    U = [USubpb1;zeros(2*size(EnrHAndTipEleIndex,1),1) ];
    W = [udual;zeros(4*size(EnrHAndTipEleIndex,1),1) ];
    v = [vdual;zeros(2*size(EnrHAndTipEleIndex,1),1) ];
end

% ====== Initialize variables ======
Uhat = U; U = [U;zeros(DIM*NodesNumPerEle,1)]; v = [v;zeros(DIM*NodesNumPerEle,1)]; 
F = [F;zeros(DIM^2*NodesNumPerEle,1)]; W = [W;zeros(DIM^2*NodesNumPerEle,1)];
% FMinusW1 = F(1:2:end)-W(1:2:end); FMinusW2 = F(2:2:end)-W(2:2:end); % Comment old codes
UMinusv = U-v; FMinusW = F-W;
CrackTipOrNot = 1-CrackTipOrNot; % Keep consistent with old version codes


% ====== Gaussian quadrature parameter ======
% ------ 3*3 Gaussian points ------
gqpt1 = 0; gqpt2 = sqrt(3/5); gqpt3 = -sqrt(3/5); gqpt = [gqpt1,gqpt2,gqpt3];
gqwt1 = 8/9; gqwt2 = 5/9; gqwt3 = 5/9; gqwt = [gqwt1,gqwt2,gqwt3];
% ------ 4*4 Gaussian points ------
% gqpt1 = 0.339981; gqpt2 = -0.339981; gqpt3 = 0.861136; gqpt4 = -0.861136;  
% gqwt1 = 0.652145; gqwt2 = 0.652145; gqwt3 = 0.347855; gqwt4 = 0.347855;  
% gqpt = [gqpt1,gqpt2,gqpt3,gqpt4]; gqwt = [gqwt1,gqwt2,gqwt3,gqwt4];
% ------ 5*5 Gaussian points ------
% gqpt1 = 0; gqpt2 = 0.538469; gqpt3 = -0.538469; gqpt4 = 0.90618; gqpt5 = -0.90618;
% gqwt1 = 0.568889; gqwt2 = 0.478629; gqwt3 = 0.478629; gqwt4 = 0.236927; gqwt5 = 0.236927;
% gqpt = [gqpt1,gqpt2,gqpt3,gqpt4,gqpt5]; gqwt = [gqwt1,gqwt2,gqwt3,gqwt4,gqwt5];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== Initialize Global FEM solver ======
ConvergeOrNot = 0; IterStep = 0;


while ConvergeOrNot < 0.5 && IterStep < 1 % Only do one step; IterStep will add value "1" in the first iteration.
    
    IterStep = IterStep+1;
    
    if clusterNo== 0 || clusterNo== 1
        hbar = waitbar(0, ['FE-method to solve Subproblem 2.']);
    else
        hbar=parfor_progressbar(size(elementsFEM,1),['FE-method to solve Subproblem 2.']);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ============= Each element, assemble stiffness matrix ============
    parfor j = 1:size(elementsFEM,1)
         
        if clusterNo== 0 || clusterNo== 1
            waitbar(j/size(elementsFEM,1));
        else
            hbar.iterate(1);
        end
        
        tempA = zeros(DIM*8,DIM*8); tempb = tempA(:,1);
        
        % ------ Find four corner points ------
        point1x = coordinatesFEM(elementsFEM(j,1),1); point1y = coordinatesFEM(elementsFEM(j,1),2);
        point2x = coordinatesFEM(elementsFEM(j,2),1); point2y = coordinatesFEM(elementsFEM(j,2),2);
        point3x = coordinatesFEM(elementsFEM(j,3),1); point3y = coordinatesFEM(elementsFEM(j,3),2);
        point4x = coordinatesFEM(elementsFEM(j,4),1); point4y = coordinatesFEM(elementsFEM(j,4),2);
        
        % ------ Find mid points 5/6/7/8 -------
        if elementsFEM(j,5) ~= 0
            point5x = coordinatesFEM(elementsFEM(j,5),1); point5y = coordinatesFEM(elementsFEM(j,5),2);
        else, point5x = 0; point5y = 0; end
        if elementsFEM(j,6) ~= 0
            point6x = coordinatesFEM(elementsFEM(j,6),1); point6y = coordinatesFEM(elementsFEM(j,6),2);
        else, point6x = 0; point6y = 0; end
        if elementsFEM(j,7) ~= 0
            point7x = coordinatesFEM(elementsFEM(j,7),1); point7y = coordinatesFEM(elementsFEM(j,7),2);
        else, point7x = 0; point7y = 0; end
        if elementsFEM(j,8) ~= 0
            point8x = coordinatesFEM(elementsFEM(j,8),1); point8y = coordinatesFEM(elementsFEM(j,8),2);
        else, point8x = 0; point8y = 0; end
        
        % ------ Calculate ksi and eta ------
        % lMatrix = [ point1x*point1y point1x point1y 1;
        %             point2x*point2y point2x point2y 1;
        %             point3x*point3y point3x point3y 1;
        %             point4x*point4y point4x point4y 1 ];
          
        % % ------ Find the element nodal indices ------
        % lb = [-1;1;1;-1]; l = linsolve(lMatrix,lb);
        % mb = [-1;-1;1;1]; m = linsolve(lMatrix,mb);
        
        % ------ Find the element nodal indices ------
        tempIndexU = [2*elementsFEM(j,1)-1 2*elementsFEM(j,1) 2*elementsFEM(j,2)-1 2*elementsFEM(j,2) ...
                2*elementsFEM(j,3)-1 2*elementsFEM(j,3) 2*elementsFEM(j,4)-1 2*elementsFEM(j,4) ...
                2*elementsFEM(j,5)-1 2*elementsFEM(j,5) 2*elementsFEM(j,6)-1 2*elementsFEM(j,6) ...
                2*elementsFEM(j,7)-1 2*elementsFEM(j,7) 2*elementsFEM(j,8)-1 2*elementsFEM(j,8)];
        
%         tp = ones(1,DIM);
%         tempIndexU = 2*elementsFEM(j,[tp,2*tp,3*tp,4*tp,5*tp,6*tp,7*tp,8*tp]);
%         tempIndexU(1:2:end) = tempIndexU(1:2:end)-1;  % size of tempIndexU: 1*16

        tempIndexF = [4*elementsFEM(j,1)-3 4*elementsFEM(j,1)-1 4*elementsFEM(j,1)-2 4*elementsFEM(j,1);
                4*elementsFEM(j,1)-3 4*elementsFEM(j,1)-1 4*elementsFEM(j,1)-2 4*elementsFEM(j,1);
                4*elementsFEM(j,2)-3 4*elementsFEM(j,2)-1 4*elementsFEM(j,2)-2 4*elementsFEM(j,2);
                4*elementsFEM(j,2)-3 4*elementsFEM(j,2)-1 4*elementsFEM(j,2)-2 4*elementsFEM(j,2);
                4*elementsFEM(j,3)-3 4*elementsFEM(j,3)-1 4*elementsFEM(j,3)-2 4*elementsFEM(j,3);
                4*elementsFEM(j,3)-3 4*elementsFEM(j,3)-1 4*elementsFEM(j,3)-2 4*elementsFEM(j,3);
                4*elementsFEM(j,4)-3 4*elementsFEM(j,4)-1 4*elementsFEM(j,4)-2 4*elementsFEM(j,4);
                4*elementsFEM(j,4)-3 4*elementsFEM(j,4)-1 4*elementsFEM(j,4)-2 4*elementsFEM(j,4);
                4*elementsFEM(j,5)-3 4*elementsFEM(j,5)-1 4*elementsFEM(j,5)-2 4*elementsFEM(j,5);
                4*elementsFEM(j,5)-3 4*elementsFEM(j,5)-1 4*elementsFEM(j,5)-2 4*elementsFEM(j,5);
                4*elementsFEM(j,6)-3 4*elementsFEM(j,6)-1 4*elementsFEM(j,6)-2 4*elementsFEM(j,6);
                4*elementsFEM(j,6)-3 4*elementsFEM(j,6)-1 4*elementsFEM(j,6)-2 4*elementsFEM(j,6);
                4*elementsFEM(j,7)-3 4*elementsFEM(j,7)-1 4*elementsFEM(j,7)-2 4*elementsFEM(j,7);
                4*elementsFEM(j,7)-3 4*elementsFEM(j,7)-1 4*elementsFEM(j,7)-2 4*elementsFEM(j,7);
                4*elementsFEM(j,8)-3 4*elementsFEM(j,8)-1 4*elementsFEM(j,8)-2 4*elementsFEM(j,8);
                4*elementsFEM(j,8)-3 4*elementsFEM(j,8)-1 4*elementsFEM(j,8)-2 4*elementsFEM(j,8)];
        
        % We don't want temp <= 0, instead, put them to the end
        for tempi = 5:8
            if elementsFEM(j,tempi) == 0
                tempIndexU(2*tempi-1)   = 2*(FEMSize/DIM+1)-1;   tempIndexU(2*tempi)   = 2*(FEMSize/DIM+1);
                tempIndexF(2*tempi-1,1) = 4*(FEMSize/DIM+1)-3;   tempIndexF(2*tempi,1) = 4*(FEMSize/DIM+1)-3;
                tempIndexF(2*tempi-1,2) = 4*(FEMSize/DIM+1)-1;   tempIndexF(2*tempi,2) = 4*(FEMSize/DIM+1)-1;
                tempIndexF(2*tempi-1,3) = 4*(FEMSize/DIM+1)-2;   tempIndexF(2*tempi,3) = 4*(FEMSize/DIM+1)-2;
                tempIndexF(2*tempi-1,4) = 4*(FEMSize/DIM+1);     tempIndexF(2*tempi,4) = 4*(FEMSize/DIM+1);
            end
        end
        
        % ------ Find the enriched functions nodal indices -------
%         if CrackOrNot == 1  % && elementsFEM(j,end) ~= 0 % Crack with and without crack tip
%             NHIndexNode = (FEMSize/DIM+1)*ones(8,1); NTipIndexNode = (FEMSize/DIM+1)*ones(8,1);
%             for tempi = 1:8
%                 if elementsFEM(j,tempi) ~= 0
%                     [tempindex,~] = find(EnrHAndTipEleIndex(:,1)==elementsFEM(j,tempi));
%                     if isempty(tempindex) ~= 1
%                         NHIndexNode(tempi) = EnrHAndTipEleIndex(tempindex,2);
%                         EleCrackHOrNot = EleCrackHOrNot + 1;
%                     end
%                     if CrackTipOrNot == 0
%                     [tempindexTip,~] = find(EnrTipEleIndex(:,1)==elementsFEM(j,tempi));
%                         if isempty(tempindexTip) ~= 1
%                             NTipIndexNode(tempi) = EnrTipEleIndex(tempindexTip,2);
%                             EleCrackTipOrNot = EleCrackTipOrNot + 1;
%                         end
%                     end
%                 end
%             end
%             
%             if EleCrackHOrNot > 0 % == 4+deltaPt5+deltaPt6+deltaPt7+deltaPt8
%                 tempNH = [  2*NHIndexNode(1)-1 2*NHIndexNode(1) 2*NHIndexNode(2)-1 2*NHIndexNode(2) ...
%                     2*NHIndexNode(3)-1 2*NHIndexNode(3) 2*NHIndexNode(4)-1 2*NHIndexNode(4) ...
%                     2*NHIndexNode(5)-1 2*NHIndexNode(5) 2*NHIndexNode(6)-1 2*NHIndexNode(6) ...
%                     2*NHIndexNode(7)-1 2*NHIndexNode(7) 2*NHIndexNode(8)-1 2*NHIndexNode(8)];
%                 
%                 tempNHF = [ 4*NHIndexNode(1)-3 4*NHIndexNode(1)-1 4*NHIndexNode(1)-2 4*NHIndexNode(1);
%                     4*NHIndexNode(1)-3 4*NHIndexNode(1)-1 4*NHIndexNode(1)-2 4*NHIndexNode(1);
%                     4*NHIndexNode(2)-3 4*NHIndexNode(2)-1 4*NHIndexNode(2)-2 4*NHIndexNode(2);
%                     4*NHIndexNode(2)-3 4*NHIndexNode(2)-1 4*NHIndexNode(2)-2 4*NHIndexNode(2);
%                     4*NHIndexNode(3)-3 4*NHIndexNode(3)-1 4*NHIndexNode(3)-2 4*NHIndexNode(3);
%                     4*NHIndexNode(3)-3 4*NHIndexNode(3)-1 4*NHIndexNode(3)-2 4*NHIndexNode(3);
%                     4*NHIndexNode(4)-3 4*NHIndexNode(4)-1 4*NHIndexNode(4)-2 4*NHIndexNode(4);
%                     4*NHIndexNode(4)-3 4*NHIndexNode(4)-1 4*NHIndexNode(4)-2 4*NHIndexNode(4);
%                     4*NHIndexNode(5)-3 4*NHIndexNode(5)-1 4*NHIndexNode(5)-2 4*NHIndexNode(5);
%                     4*NHIndexNode(5)-3 4*NHIndexNode(5)-1 4*NHIndexNode(5)-2 4*NHIndexNode(5);
%                     4*NHIndexNode(6)-3 4*NHIndexNode(6)-1 4*NHIndexNode(6)-2 4*NHIndexNode(6);
%                     4*NHIndexNode(6)-3 4*NHIndexNode(6)-1 4*NHIndexNode(6)-2 4*NHIndexNode(6);
%                     4*NHIndexNode(7)-3 4*NHIndexNode(7)-1 4*NHIndexNode(7)-2 4*NHIndexNode(7);
%                     4*NHIndexNode(7)-3 4*NHIndexNode(7)-1 4*NHIndexNode(7)-2 4*NHIndexNode(7);
%                     4*NHIndexNode(8)-3 4*NHIndexNode(8)-1 4*NHIndexNode(8)-2 4*NHIndexNode(8);
%                     4*NHIndexNode(8)-3 4*NHIndexNode(8)-1 4*NHIndexNode(8)-2 4*NHIndexNode(8)];
%                 
%                 tempIndexU = [tempIndexU tempNH]; tempIndexF = [tempIndexF;tempNHF];
%                 
%                 if EleCrackTipOrNot > 0 
%                     
%                     tempNTip = [2*NTipIndexNode(1)-1 2*NTipIndexNode(1) 2*(NTipIndexNode(1)+1)-1 2*(NTipIndexNode(1)+1) ...
%                         2*(NTipIndexNode(1)+2)-1 2*(NTipIndexNode(1)+2) 2*(NTipIndexNode(1)+3)-1 2*(NTipIndexNode(1)+3) ...
%                         2*NTipIndexNode(2)-1 2*NTipIndexNode(2) 2*(NTipIndexNode(2)+1)-1 2*(NTipIndexNode(2)+1) ...
%                         2*(NTipIndexNode(2)+2)-1 2*(NTipIndexNode(2)+2) 2*(NTipIndexNode(2)+3)-1 2*(NTipIndexNode(2)+3) ...
%                         2*NTipIndexNode(3)-1 2*NTipIndexNode(3) 2*(NTipIndexNode(3)+1)-1 2*(NTipIndexNode(3)+1) ...
%                         2*(NTipIndexNode(3)+2)-1 2*(NTipIndexNode(3)+2) 2*(NTipIndexNode(3)+3)-1 2*(NTipIndexNode(3)+3) ...
%                         2*NTipIndexNode(4)-1 2*NTipIndexNode(4) 2*(NTipIndexNode(4)+1)-1 2*(NTipIndexNode(4)+1) ...
%                         2*(NTipIndexNode(4)+2)-1 2*(NTipIndexNode(4)+2) 2*(NTipIndexNode(4)+3)-1 2*(NTipIndexNode(4)+3) ...
%                         2*NTipIndexNode(5)-1 2*NTipIndexNode(5) 2*(NTipIndexNode(5)+1)-1 2*(NTipIndexNode(5)+1) ...
%                         2*(NTipIndexNode(5)+2)-1 2*(NTipIndexNode(5)+2) 2*(NTipIndexNode(5)+3)-1 2*(NTipIndexNode(5)+3) ...
%                         2*NTipIndexNode(6)-1 2*NTipIndexNode(6)  2*(NTipIndexNode(6)+1)-1 2*(NTipIndexNode(6)+1) ...
%                         2*(NTipIndexNode(6)+2)-1 2*(NTipIndexNode(6)+2)  2*(NTipIndexNode(6)+3)-1 2*(NTipIndexNode(6)+3) ...
%                         2*NTipIndexNode(7)-1 2*NTipIndexNode(7) 2*(NTipIndexNode(7)+1)-1 2*(NTipIndexNode(7)+1) ...
%                         2*(NTipIndexNode(7)+2)-1 2*(NTipIndexNode(7)+2)  2*(NTipIndexNode(7)+3)-1 2*(NTipIndexNode(7)+3) ...
%                         2*NTipIndexNode(8)-1 2*NTipIndexNode(8) 2*(NTipIndexNode(8)+1)-1 2*(NTipIndexNode(8)+1) ...
%                         2*(NTipIndexNode(8)+2)-1 2*(NTipIndexNode(8)+2) 2*(NTipIndexNode(8)+3)-1 2*(NTipIndexNode(8)+3) ];
%                      
%                     tempNTipF = [4*NTipIndexNode(1)-3 4*NTipIndexNode(1)-1 4*NTipIndexNode(1)-2 4*NTipIndexNode(1);
%                         4*NTipIndexNode(1)-3 4*NTipIndexNode(1)-1 4*NTipIndexNode(1)-2 4*NTipIndexNode(1);
%                         4*(NTipIndexNode(1)+1)-3 4*(NTipIndexNode(1)+1)-1 4*(NTipIndexNode(1)+1)-2 4*(NTipIndexNode(1)+1);
%                         4*(NTipIndexNode(1)+1)-3 4*(NTipIndexNode(1)+1)-1 4*(NTipIndexNode(1)+1)-2 4*(NTipIndexNode(1)+1);
%                         4*(NTipIndexNode(1)+2)-3 4*(NTipIndexNode(1)+2)-1 4*(NTipIndexNode(1)+2)-2 4*(NTipIndexNode(1)+2);
%                         4*(NTipIndexNode(1)+2)-3 4*(NTipIndexNode(1)+2)-1 4*(NTipIndexNode(1)+2)-2 4*(NTipIndexNode(1)+2);
%                         4*(NTipIndexNode(1)+3)-3 4*(NTipIndexNode(1)+3)-1 4*(NTipIndexNode(1)+3)-2 4*(NTipIndexNode(1)+3);
%                         4*(NTipIndexNode(1)+3)-3 4*(NTipIndexNode(1)+3)-1 4*(NTipIndexNode(1)+3)-2 4*(NTipIndexNode(1)+3);
%                         4*NTipIndexNode(2)-3 4*NTipIndexNode(2)-1 4*NTipIndexNode(2)-2 4*NTipIndexNode(2);
%                         4*NTipIndexNode(2)-3 4*NTipIndexNode(2)-1 4*NTipIndexNode(2)-2 4*NTipIndexNode(2);
%                         4*(NTipIndexNode(2)+1)-3 4*(NTipIndexNode(2)+1)-1 4*(NTipIndexNode(2)+1)-2 4*(NTipIndexNode(2)+1);
%                         4*(NTipIndexNode(2)+1)-3 4*(NTipIndexNode(2)+1)-1 4*(NTipIndexNode(2)+1)-2 4*(NTipIndexNode(2)+1);
%                         4*(NTipIndexNode(2)+2)-3 4*(NTipIndexNode(2)+2)-1 4*(NTipIndexNode(2)+2)-2 4*(NTipIndexNode(2)+2);
%                         4*(NTipIndexNode(2)+2)-3 4*(NTipIndexNode(2)+2)-1 4*(NTipIndexNode(2)+2)-2 4*(NTipIndexNode(2)+2);
%                         4*(NTipIndexNode(2)+3)-3 4*(NTipIndexNode(2)+3)-1 4*(NTipIndexNode(2)+3)-2 4*(NTipIndexNode(2)+3);
%                         4*(NTipIndexNode(2)+3)-3 4*(NTipIndexNode(2)+3)-1 4*(NTipIndexNode(2)+3)-2 4*(NTipIndexNode(2)+3);
%                         4*NTipIndexNode(3)-3 4*NTipIndexNode(3)-1 4*NTipIndexNode(3)-2 4*NTipIndexNode(3);
%                         4*NTipIndexNode(3)-3 4*NTipIndexNode(3)-1 4*NTipIndexNode(3)-2 4*NTipIndexNode(3);
%                         4*(NTipIndexNode(3)+1)-3 4*(NTipIndexNode(3)+1)-1 4*(NTipIndexNode(3)+1)-2 4*(NTipIndexNode(3)+1);
%                         4*(NTipIndexNode(3)+1)-3 4*(NTipIndexNode(3)+1)-1 4*(NTipIndexNode(3)+1)-2 4*(NTipIndexNode(3)+1);
%                         4*(NTipIndexNode(3)+2)-3 4*(NTipIndexNode(3)+2)-1 4*(NTipIndexNode(3)+2)-2 4*(NTipIndexNode(3)+2);
%                         4*(NTipIndexNode(3)+2)-3 4*(NTipIndexNode(3)+2)-1 4*(NTipIndexNode(3)+2)-2 4*(NTipIndexNode(3)+2);
%                         4*(NTipIndexNode(3)+3)-3 4*(NTipIndexNode(3)+3)-1 4*(NTipIndexNode(3)+3)-2 4*(NTipIndexNode(3)+3);
%                         4*(NTipIndexNode(3)+3)-3 4*(NTipIndexNode(3)+3)-1 4*(NTipIndexNode(3)+3)-2 4*(NTipIndexNode(3)+3);
%                         4*NTipIndexNode(4)-3 4*NTipIndexNode(4)-1 4*NTipIndexNode(4)-2 4*NTipIndexNode(4);
%                         4*NTipIndexNode(4)-3 4*NTipIndexNode(4)-1 4*NTipIndexNode(4)-2 4*NTipIndexNode(4);
%                         4*(NTipIndexNode(4)+1)-3 4*(NTipIndexNode(4)+1)-1 4*(NTipIndexNode(4)+1)-2 4*(NTipIndexNode(4)+1);
%                         4*(NTipIndexNode(4)+1)-3 4*(NTipIndexNode(4)+1)-1 4*(NTipIndexNode(4)+1)-2 4*(NTipIndexNode(4)+1);
%                         4*(NTipIndexNode(4)+2)-3 4*(NTipIndexNode(4)+2)-1 4*(NTipIndexNode(4)+2)-2 4*(NTipIndexNode(4)+2);
%                         4*(NTipIndexNode(4)+2)-3 4*(NTipIndexNode(4)+2)-1 4*(NTipIndexNode(4)+2)-2 4*(NTipIndexNode(4)+2);
%                         4*(NTipIndexNode(4)+3)-3 4*(NTipIndexNode(4)+3)-1 4*(NTipIndexNode(4)+3)-2 4*(NTipIndexNode(4)+3);
%                         4*(NTipIndexNode(4)+3)-3 4*(NTipIndexNode(4)+3)-1 4*(NTipIndexNode(4)+3)-2 4*(NTipIndexNode(4)+3);
%                         4*NTipIndexNode(5)-3 4*NTipIndexNode(5)-1 4*NTipIndexNode(5)-2 4*NTipIndexNode(5);
%                         4*NTipIndexNode(5)-3 4*NTipIndexNode(5)-1 4*NTipIndexNode(5)-2 4*NTipIndexNode(5);
%                         4*(NTipIndexNode(5)+1)-3 4*(NTipIndexNode(5)+1)-1 4*(NTipIndexNode(5)+1)-2 4*(NTipIndexNode(5)+1);
%                         4*(NTipIndexNode(5)+1)-3 4*(NTipIndexNode(5)+1)-1 4*(NTipIndexNode(5)+1)-2 4*(NTipIndexNode(5)+1);
%                         4*(NTipIndexNode(5)+2)-3 4*(NTipIndexNode(5)+2)-1 4*(NTipIndexNode(5)+2)-2 4*(NTipIndexNode(5)+2);
%                         4*(NTipIndexNode(5)+2)-3 4*(NTipIndexNode(5)+2)-1 4*(NTipIndexNode(5)+2)-2 4*(NTipIndexNode(5)+2);
%                         4*(NTipIndexNode(5)+3)-3 4*(NTipIndexNode(5)+3)-1 4*(NTipIndexNode(5)+3)-2 4*(NTipIndexNode(5)+3);
%                         4*(NTipIndexNode(5)+3)-3 4*(NTipIndexNode(5)+3)-1 4*(NTipIndexNode(5)+3)-2 4*(NTipIndexNode(5)+3);
%                         4*NTipIndexNode(6)-3 4*NTipIndexNode(6)-1 4*NTipIndexNode(6)-2 4*NTipIndexNode(6);
%                         4*NTipIndexNode(6)-3 4*NTipIndexNode(6)-1 4*NTipIndexNode(6)-2 4*NTipIndexNode(6);
%                         4*(NTipIndexNode(6)+1)-3 4*(NTipIndexNode(6)+1)-1 4*(NTipIndexNode(6)+1)-2 4*(NTipIndexNode(6)+1);
%                         4*(NTipIndexNode(6)+1)-3 4*(NTipIndexNode(6)+1)-1 4*(NTipIndexNode(6)+1)-2 4*(NTipIndexNode(6)+1);
%                         4*(NTipIndexNode(6)+2)-3 4*(NTipIndexNode(6)+2)-1 4*(NTipIndexNode(6)+2)-2 4*(NTipIndexNode(6)+2);
%                         4*(NTipIndexNode(6)+2)-3 4*(NTipIndexNode(6)+2)-1 4*(NTipIndexNode(6)+2)-2 4*(NTipIndexNode(6)+2);
%                         4*(NTipIndexNode(6)+3)-3 4*(NTipIndexNode(6)+3)-1 4*(NTipIndexNode(6)+3)-2 4*(NTipIndexNode(6)+3);
%                         4*(NTipIndexNode(6)+3)-3 4*(NTipIndexNode(6)+3)-1 4*(NTipIndexNode(6)+3)-2 4*(NTipIndexNode(6)+3);
%                         4*NTipIndexNode(7)-3 4*NTipIndexNode(7)-1 4*NTipIndexNode(7)-2 4*NTipIndexNode(7);
%                         4*NTipIndexNode(7)-3 4*NTipIndexNode(7)-1 4*NTipIndexNode(7)-2 4*NTipIndexNode(7);
%                         4*(NTipIndexNode(7)+1)-3 4*(NTipIndexNode(7)+1)-1 4*(NTipIndexNode(7)+1)-2 4*(NTipIndexNode(7)+1);
%                         4*(NTipIndexNode(7)+1)-3 4*(NTipIndexNode(7)+1)-1 4*(NTipIndexNode(7)+1)-2 4*(NTipIndexNode(7)+1);
%                         4*(NTipIndexNode(7)+2)-3 4*(NTipIndexNode(7)+2)-1 4*(NTipIndexNode(7)+2)-2 4*(NTipIndexNode(7)+2);
%                         4*(NTipIndexNode(7)+2)-3 4*(NTipIndexNode(7)+2)-1 4*(NTipIndexNode(7)+2)-2 4*(NTipIndexNode(7)+2);
%                         4*(NTipIndexNode(7)+3)-3 4*(NTipIndexNode(7)+3)-1 4*(NTipIndexNode(7)+3)-2 4*(NTipIndexNode(7)+3);
%                         4*(NTipIndexNode(7)+3)-3 4*(NTipIndexNode(7)+3)-1 4*(NTipIndexNode(7)+3)-2 4*(NTipIndexNode(7)+3);
%                         4*NTipIndexNode(8)-3 4*NTipIndexNode(8)-1 4*NTipIndexNode(8)-2 4*NTipIndexNode(8);
%                         4*NTipIndexNode(8)-3 4*NTipIndexNode(8)-1 4*NTipIndexNode(8)-2 4*NTipIndexNode(8);
%                         4*(NTipIndexNode(8)+1)-3 4*(NTipIndexNode(8)+1)-1 4*(NTipIndexNode(8)+1)-2 4*(NTipIndexNode(8)+1);
%                         4*(NTipIndexNode(8)+1)-3 4*(NTipIndexNode(8)+1)-1 4*(NTipIndexNode(8)+1)-2 4*(NTipIndexNode(8)+1);
%                         4*(NTipIndexNode(8)+2)-3 4*(NTipIndexNode(8)+2)-1 4*(NTipIndexNode(8)+2)-2 4*(NTipIndexNode(8)+2);
%                         4*(NTipIndexNode(8)+2)-3 4*(NTipIndexNode(8)+2)-1 4*(NTipIndexNode(8)+2)-2 4*(NTipIndexNode(8)+2);
%                         4*(NTipIndexNode(8)+3)-3 4*(NTipIndexNode(8)+3)-1 4*(NTipIndexNode(8)+3)-2 4*(NTipIndexNode(8)+3);
%                         4*(NTipIndexNode(8)+3)-3 4*(NTipIndexNode(8)+3)-1 4*(NTipIndexNode(8)+3)-2 4*(NTipIndexNode(8)+3)];
%                     
%                     tempIndexU = [tempIndexU tempNTip]; tempIndexF = [tempIndexF;tempNTipF];
%                      
%                 end
%             end
%             
%         end
        
        % ------ Set Gauss points ------
        % pt1 = elementsFEM(j,1); pt2 = elementsFEM(j,3);
        % pointOfx = zeros(length(gqwt),1); pointOfy = zeros(length(gqwt),1);
        % for tempi = 1:length(gqwt)
        %     pointOfx(tempi) = gqpt(tempi)*0.5*(coordinatesFEM(pt2,1)-coordinatesFEM(pt1,1))+0.5*(coordinatesFEM(pt2,1)+coordinatesFEM(pt1,1));
        %     pointOfy(tempi) = gqpt(tempi)*0.5*(coordinatesFEM(pt2,2)-coordinatesFEM(pt1,2))+0.5*(coordinatesFEM(pt2,2)+coordinatesFEM(pt1,2));
        % end
        %
        % for tempi = 1:size(pointOfx,1)
        %     for tempj = 1:size(pointOfy,1)
        % % ------ Calculate ksi and eta ------
        %        ksi = l(1)*pointOfx(tempi)*pointOfy(tempj) + l(2)*pointOfx(tempi) + l(3)*pointOfy(tempj) + l(4) ;
        %        eta = m(1)*pointOfx(tempi)*pointOfy(tempj) + m(2)*pointOfx(tempi) + m(3)*pointOfy(tempj) + m(4) ;
        % %%%%%% Use following codes %%%%%%
        % ------ Nine Gauss integral points ------
        ksietaList = [-sqrt(3/5), 0, sqrt(3/5)]; 
        ksietaWeightList = [5/9, 8/9, 5/9];
        %ksietaList = [-0.90618,-0.538469,0,0.538469,0.90618]; 
        %ksietaWeightList = [0.236927,0.478629,0.568889,0.478629,0.236927];
        
        [ksiAll, etaAll] = ndgrid(ksietaList, ksietaList);
        [weightksiAll, weightetaAll] = ndgrid(ksietaWeightList, ksietaWeightList);
        
        ksiAll = ksiAll(:); etaAll = etaAll(:); weightksiAll = weightksiAll(:); weightetaAll = weightetaAll(:);
        
        
        %for tempi = 1:length(ksietaList)
        %    for tempj = 1:length(ksietaList)
        for tempjj = 1:length(ksiAll)
            
            % ---------------------------
               % ksi = ksietaList(tempi);
               % eta = ksietaList(tempj);
               % weightksi = ksietaWeightList(tempi);
               % weighteta = ksietaWeightList(tempj);
               ksi = ksiAll(tempjj); eta = etaAll(tempjj); weightksi = weightksiAll(tempjj); weighteta = weightetaAll(tempjj);
               
         
                % ------ Calculate N ------
                deltaPt5 = 0; deltaPt6 = 0; deltaPt7 = 0; deltaPt8 = 0;
                if elementsFEM(j,5) ~= 0; deltaPt5 = 1; end
                if elementsFEM(j,6) ~= 0; deltaPt6 = 1; end
                if elementsFEM(j,7) ~= 0; deltaPt7 = 1; end
                if elementsFEM(j,8) ~= 0; deltaPt8 = 1; end
               
                N5 = deltaPt5*0.5*(1+ksi)*(1-eta^2);
                N6 = deltaPt6*0.5*(1+eta)*(1-ksi^2);
                N7 = deltaPt7*0.5*(1-ksi)*(1-eta^2);
                N8 = deltaPt8*0.5*(1-eta)*(1-ksi^2);
                %N5 = deltaPt5*0.5*(1+ksi)*(1-abs(eta));
                %N6 = deltaPt6*0.5*(1+eta)*(1-abs(ksi));
                %N7 = deltaPt7*0.5*(1-ksi)*(1-abs(eta));
                %N8 = deltaPt8*0.5*(1-eta)*(1-abs(ksi));
                
                N1 = (1-ksi)*(1-eta)*0.25 - 0.5*(N7+N8);
                N2 = (1+ksi)*(1-eta)*0.25 - 0.5*(N8+N5);
                N3 = (1+ksi)*(1+eta)*0.25 - 0.5*(N5+N6);
                N4 = (1-ksi)*(1+eta)*0.25 - 0.5*(N6+N7);
                
                %%%%%%%%% Consider crack enriched basis %%%%%%%%%%%%%%%
                if CrackOrNot == 0 || EleCrackHOrNot == 0 % No crack
                    H1 = 0; H2 = 0; H3 = 0; H4 = 0; H5 = 0; H6 = 0; H7 = 0; H8 = 0; 
                    C1 = 0; C2 = 0; C3 = 0; C4 = 0;
                elseif CrackOrNot == 1 && EleCrackHOrNot > 0 % == 4+deltaPt5+deltaPt6+deltaPt7+deltaPt8 % Crack both with and without crack tip
                    % Compute NH function: NH = N*H; so we need to compute
                    % H function, to compute H function, we need to compare
                    % coordinatesFEM of Gauss points with other nodes
                    temp1 = CrackPathCen(1)*(point1x) + CrackPathCen(2)*(point1y) + 1; % Node 1
                    temp2 = CrackPathCen(1)*(point2x) + CrackPathCen(2)*(point2y) + 1; % Node 2
                    temp3 = CrackPathCen(1)*(point3x) + CrackPathCen(2)*(point3y) + 1; % Node 3
                    temp4 = CrackPathCen(1)*(point4x) + CrackPathCen(2)*(point4y) + 1; % Node 4
                    temp5 = deltaPt5 * (CrackPathCen(1)*(point5x) + CrackPathCen(2)*(point5y) + 1); % Node 5
                    temp6 = deltaPt6 * (CrackPathCen(1)*(point6x) + CrackPathCen(2)*(point6y) + 1); % Node 6
                    temp7 = deltaPt7 * (CrackPathCen(1)*(point7x) + CrackPathCen(2)*(point7y) + 1); % Node 7
                    temp8 = deltaPt8 * (CrackPathCen(1)*(point8x) + CrackPathCen(2)*(point8y) + 1); % Node 8
                    tempNH = CrackPathCen(1)*pointOfx(tempi) + CrackPathCen(2)*pointOfy(tempj) + 1;
                    % If temp * tempNH < 0, then they are not at the same sides
                    H1 = 0; if (temp1 * tempNH < 0); H1 = sign(temp1-tempNH); end
                    H2 = 0; if (temp2 * tempNH < 0); H2 = sign(temp2-tempNH); end 
                    H3 = 0; if (temp3 * tempNH < 0); H3 = sign(temp3-tempNH); end
                    H4 = 0; if (temp4 * tempNH < 0); H4 = sign(temp4-tempNH); end                     
                    H5 = 0; if (temp5 * tempNH < 0); H5 = sign(temp5-tempNH); end
                    H6 = 0; if (temp6 * tempNH < 0); H6 = sign(temp6-tempNH); end                    
                    H7 = 0; if (temp7 * tempNH < 0); H7 = sign(temp7-tempNH); end
                    H8 = 0; if (temp8 * tempNH < 0); H8 = sign(temp8-tempNH); end  
                    % [tempNH, temp1, temp2, temp3, temp4, CrackTip(1), pointOfx(tempi), H1, H2, H3, H4 ]
                    
                    C1 = 0; C2 = 0; C3 = 0; C4 = 0; 
                    if EleCrackTipOrNot > 0 % if elementsFEM(j,end) == 2 % Crack with crack tip
                        % Compute C11-C44 function: Cij = Ni*Cj;
                        tempr = (sqrt( (pointOfx(tempi)-CrackTip(1))^2 + (pointOfy(tempj)-CrackTip(2))^2 ));
                        temptheta = atan2( (pointOfx(tempi)-CrackTip(1)) , (pointOfy(tempj)-CrackTip(2)) );  % atan2(Y,X)
                        C1 = sqrt(tempr)*sin(0.5*temptheta); 
                        C2 = sqrt(tempr)*cos(0.5*temptheta);
                        C3 = sqrt(tempr)*sin(0.5*temptheta)*sin(temptheta); 
                        C4 = sqrt(tempr)*cos(0.5*temptheta)*sin(temptheta);
                         
                        deltaCNode1 = 0; deltaCNode2 = 0; deltaCNode3 = 0; deltaCNode4 = 0; deltaCNode5 = 0; deltaCNode6 = 0; deltaCNode7 = 0; deltaCNode8 = 0;
                        if NTipIndexNode(1) ~= (FEMSize/DIM+1); deltaCNode1 = 1; end
                        if NTipIndexNode(2) ~= (FEMSize/DIM+1); deltaCNode2 = 1; end
                        if NTipIndexNode(3) ~= (FEMSize/DIM+1); deltaCNode3 = 1; end
                        if NTipIndexNode(4) ~= (FEMSize/DIM+1); deltaCNode4 = 1; end
                        if NTipIndexNode(5) ~= (FEMSize/DIM+1); deltaCNode5 = 1; end
                        if NTipIndexNode(6) ~= (FEMSize/DIM+1); deltaCNode6 = 1; end
                        if NTipIndexNode(7) ~= (FEMSize/DIM+1); deltaCNode7 = 1; end
                        if NTipIndexNode(8) ~= (FEMSize/DIM+1); deltaCNode8 = 1; end
                    
                    end
                end
                
                if CrackOrNot == 0 || EleCrackHOrNot == 0 % < 4+deltaPt5+deltaPt6+deltaPt7+deltaPt8 % No crack
                    N = [N1 0 N2 0 N3 0 N4 0 N5 0 N6 0 N7 0 N8 0;
                         0 N1 0 N2 0 N3 0 N4 0 N5 0 N6 0 N7 0 N8];
                    NDiag = diag([N1,N1,N2,N2,N3,N3,N4,N4,N5,N5,N6,N6,N7,N7,N8,N8]);
                elseif CrackOrNot == 1 && EleCrackHOrNot > 0 && EleCrackTipOrNot == 0 % == 4+deltaPt5+deltaPt6+deltaPt7+deltaPt8 && EleCrackTipOrNot == 0 % Crack but without crack tip
                    N = [N1 0 N2 0 N3 0 N4 0 N5 0 N6 0 N7 0 N8 0  ...
                         N1*H1 0 N2*H2 0 N3*H3 0 N4*H4 0 N5*H5 0 N6*H6 0 N7*H7 0 N8*H8 0;
                         0 N1 0 N2 0 N3 0 N4 0 N5 0 N6 0 N7 0 N8  ...
                         0 N1*H1 0 N2*H2 0 N3*H3 0 N4*H4 0 N5*H5 0 N6*H6 0 N7*H7 0 N8*H8];
                    NDiag = diag([N1,N1,N2,N2,N3,N3,N4,N4,N5,N5,N6,N6,N7,N7,N8,N8,...
                        N1*H1,N1*H1,N2*H2,N2*H2,N3*H3,N3*H3,N4*H4,N4*H4,N5*H5,N5*H5,N6*H6,N6*H6,N7*H7,N7*H7,N8*H8,N8*H8]);
                elseif CrackOrNot == 1 && EleCrackHOrNot > 0 && EleCrackTipOrNot > 0 % if elementsFEM(j,end) == 2 % Crack with crack tip % Crack with crack tip
                    N = [N1 0 N2 0 N3 0 N4 0 N5 0 N6 0 N7 0 N8 0  ...
                        N1*H1 0 N2*H2 0 N3*H3 0 N4*H4 0 N5*H5 0 N6*H6 0 N7*H7 0 N8*H8 0 ...
                        C1*N1*deltaCNode1 0 C2*N1*deltaCNode1 0 C3*N1*deltaCNode1 0 C4*N1*deltaCNode1 0 ...
                        C1*N2*deltaCNode2 0 C2*N2*deltaCNode2 0 C3*N2*deltaCNode2 0 C4*N2*deltaCNode2 0 ...
                        C1*N3*deltaCNode3 0 C2*N3*deltaCNode3 0 C3*N3*deltaCNode3 0 C4*N3*deltaCNode3 0 ...
                        C1*N4*deltaCNode4 0 C2*N4*deltaCNode4 0 C3*N4*deltaCNode4 0 C4*N4*deltaCNode4 0 ...
                        C1*N5*deltaCNode5 0 C2*N5*deltaCNode5 0 C3*N5*deltaCNode5 0 C4*N5*deltaCNode5 0 ...
                        C1*N6*deltaCNode6 0 C2*N6*deltaCNode6 0 C3*N6*deltaCNode6 0 C4*N6*deltaCNode6 0 ...
                        C1*N7*deltaCNode7 0 C2*N7*deltaCNode7 0 C3*N7*deltaCNode7 0 C4*N7*deltaCNode7 0 ...
                        C1*N8*deltaCNode8 0 C2*N8*deltaCNode8 0 C3*N8*deltaCNode8 0 C4*N8*deltaCNode8 0 ;
                        0 N1 0 N2 0 N3 0 N4 0 N5 0 N6 0 N7 0 N8  ...
                        0 N1*H1 0 N2*H2 0 N3*H3 0 N4*H4 0 N5*H5 0 N6*H6 0 N7*H7 0 N8*H8 ...
                        0 C1*N1*deltaCNode1 0 C2*N1*deltaCNode1 0 C3*N1*deltaCNode1 0 C4*N1*deltaCNode1 ...
                        0 C1*N2*deltaCNode2 0 C2*N2*deltaCNode2 0 C3*N2*deltaCNode2 0 C4*N2*deltaCNode2 ...
                        0 C1*N3*deltaCNode3 0 C2*N3*deltaCNode3 0 C3*N3*deltaCNode3 0 C4*N3*deltaCNode3 ...
                        0 C1*N4*deltaCNode4 0 C2*N4*deltaCNode4 0 C3*N4*deltaCNode4 0 C4*N4*deltaCNode4 ...
                        0 C1*N5*deltaCNode5 0 C2*N5*deltaCNode5 0 C3*N5*deltaCNode5 0 C4*N5*deltaCNode5 ...
                        0 C1*N6*deltaCNode6 0 C2*N6*deltaCNode6 0 C3*N6*deltaCNode6 0 C4*N6*deltaCNode6 ...
                        0 C1*N7*deltaCNode7 0 C2*N7*deltaCNode7 0 C3*N7*deltaCNode7 0 C4*N7*deltaCNode7 ...
                        0 C1*N8*deltaCNode8 0 C2*N8*deltaCNode8 0 C3*N8*deltaCNode8 0 C4*N8*deltaCNode8 ];
                    NDiag = diag([N1,N1,N2,N2,N3,N3,N4,N4,N5,N5,N6,N6,N7,N7,N8,N8,...
                        N1*H1,N1*H1,N2*H2,N2*H2,N3*H3,N3*H3,N4*H4,N4*H4,N5*H5,N5*H5,N6*H6,N6*H6,N7*H7,N7*H7,N8*H8,N8*H8,...
                        C1*N1*deltaCNode1,C1*N1*deltaCNode1,C2*N1*deltaCNode1,C2*N1*deltaCNode1,...
                        C3*N1*deltaCNode1,C3*N1*deltaCNode1,C4*N1*deltaCNode1,C4*N1*deltaCNode1,...
                        C1*N2*deltaCNode2,C1*N2*deltaCNode2,C2*N2*deltaCNode2,C2*N2*deltaCNode2,...
                        C3*N2*deltaCNode2,C3*N2*deltaCNode2,C4*N2*deltaCNode2,C4*N2*deltaCNode2,...
                        C1*N3*deltaCNode3,C1*N3*deltaCNode3,C2*N3*deltaCNode3,C2*N3*deltaCNode3,...
                        C3*N3*deltaCNode3,C3*N3*deltaCNode3,C4*N3*deltaCNode3,C4*N3*deltaCNode3,...
                        C1*N4*deltaCNode4,C1*N4*deltaCNode4,C2*N4*deltaCNode4,C2*N4*deltaCNode4,...
                        C3*N4*deltaCNode4,C3*N4*deltaCNode4,C4*N4*deltaCNode4,C4*N4*deltaCNode4,...
                        C1*N5*deltaCNode5,C1*N5*deltaCNode5,C2*N5*deltaCNode5,C2*N5*deltaCNode5,...
                        C3*N5*deltaCNode5,C3*N5*deltaCNode5,C4*N5*deltaCNode5,C4*N5*deltaCNode5,...
                        C1*N6*deltaCNode6,C1*N6*deltaCNode6,C2*N6*deltaCNode6,C2*N6*deltaCNode6,...
                        C3*N6*deltaCNode6,C3*N6*deltaCNode6,C4*N6*deltaCNode6,C4*N6*deltaCNode6,...
                        C1*N7*deltaCNode7,C1*N7*deltaCNode7,C2*N7*deltaCNode7,C2*N7*deltaCNode7,...
                        C3*N7*deltaCNode7,C3*N7*deltaCNode7,C4*N7*deltaCNode7,C4*N7*deltaCNode7,...
                        C1*N8*deltaCNode8,C1*N8*deltaCNode8,C2*N8*deltaCNode8,C2*N8*deltaCNode8,...
                        C3*N8*deltaCNode8,C3*N8*deltaCNode8,C4*N8*deltaCNode8,C4*N8*deltaCNode8]);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % NDiag = diag([N1,N1,N2,N2,N3,N3,N4,N4,N5,N5,N6,N6,N7,N7,N8,N8]);
                % ------ Build J matrix ------
                % Comment: I didn't change Jacobian matrix J when enriched 
                % functions are added.
                J11 = funDN1Dksi(ksi,eta,deltaPt7,deltaPt8)*point1x + funDN2Dksi(ksi,eta,deltaPt8,deltaPt5)*point2x + ...
                      funDN3Dksi(ksi,eta,deltaPt5,deltaPt6)*point3x + funDN4Dksi(ksi,eta,deltaPt6,deltaPt7)*point4x + ...
                      funDN5Dksi(ksi,eta,deltaPt5)*point5x + funDN6Dksi(ksi,eta,deltaPt6)*point6x + ...
                      funDN7Dksi(ksi,eta,deltaPt7)*point7x + funDN8Dksi(ksi,eta,deltaPt8)*point8x;
                J12 = funDN1Dksi(ksi,eta,deltaPt7,deltaPt8)*point1y + funDN2Dksi(ksi,eta,deltaPt8,deltaPt5)*point2y + ...
                      funDN3Dksi(ksi,eta,deltaPt5,deltaPt6)*point3y + funDN4Dksi(ksi,eta,deltaPt6,deltaPt7)*point4y + ...
                      funDN5Dksi(ksi,eta,deltaPt5)*point5y + funDN6Dksi(ksi,eta,deltaPt6)*point6y + ...
                      funDN7Dksi(ksi,eta,deltaPt7)*point7y + funDN8Dksi(ksi,eta,deltaPt8)*point8y;
                J21 = funDN1Deta(ksi,eta,deltaPt7,deltaPt8)*point1x + funDN2Deta(ksi,eta,deltaPt8,deltaPt5)*point2x + ...
                      funDN3Deta(ksi,eta,deltaPt5,deltaPt6)*point3x + funDN4Deta(ksi,eta,deltaPt6,deltaPt7)*point4x + ...
                      funDN5Deta(ksi,eta,deltaPt5)*point5x + funDN6Deta(ksi,eta,deltaPt6)*point6x + ...
                      funDN7Deta(ksi,eta,deltaPt7)*point7x + funDN8Deta(ksi,eta,deltaPt8)*point8x;
                J22 = funDN1Deta(ksi,eta,deltaPt7,deltaPt8)*point1y + funDN2Deta(ksi,eta,deltaPt8,deltaPt5)*point2y + ...
                      funDN3Deta(ksi,eta,deltaPt5,deltaPt6)*point3y + funDN4Deta(ksi,eta,deltaPt6,deltaPt7)*point4y + ...
                      funDN5Deta(ksi,eta,deltaPt5)*point5y + funDN6Deta(ksi,eta,deltaPt6)*point6y + ...
                      funDN7Deta(ksi,eta,deltaPt7)*point7y + funDN8Deta(ksi,eta,deltaPt8)*point8y;
                J = [J11 J12; J21 J22];
                
                Jacobian = det(J);
                InvJ = 1/Jacobian*[J22 -J12; -J21 J11];
                
                % ------ Compute DN matrix ------
                if CrackOrNot == 0 || EleCrackHOrNot == 0 % no crack pass
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
                    
                elseif CrackOrNot == 1 && EleCrackHOrNot > 0 && EleCrackTipOrNot == 0 % Crack but without crack tip
                    DN = [InvJ zeros(2,2); zeros(2,2) InvJ] * ...
                        [funDN1Dksi(ksi,eta,deltaPt7,deltaPt8) 0 funDN2Dksi(ksi,eta,deltaPt8,deltaPt5) 0 ...
                        funDN3Dksi(ksi,eta,deltaPt5,deltaPt6) 0 funDN4Dksi(ksi,eta,deltaPt6,deltaPt7) 0 ...
                        funDN5Dksi(ksi,eta,deltaPt5) 0 funDN6Dksi(ksi,eta,deltaPt6) 0 ...
                        funDN7Dksi(ksi,eta,deltaPt7) 0 funDN8Dksi(ksi,eta,deltaPt8) 0 ...
                        H1*funDN1Dksi(ksi,eta,deltaPt7,deltaPt8) 0 H2*funDN2Dksi(ksi,eta,deltaPt8,deltaPt5) 0 ...
                        H3*funDN3Dksi(ksi,eta,deltaPt5,deltaPt6) 0 H4*funDN4Dksi(ksi,eta,deltaPt6,deltaPt7) 0 ...
                        H5*funDN5Dksi(ksi,eta,deltaPt5) 0 H6*funDN6Dksi(ksi,eta,deltaPt6) 0 ...
                        H7*funDN7Dksi(ksi,eta,deltaPt7) 0 H8*funDN8Dksi(ksi,eta,deltaPt8) 0;
                        funDN1Deta(ksi,eta,deltaPt7,deltaPt8) 0 funDN2Deta(ksi,eta,deltaPt8,deltaPt5) 0 ...
                        funDN3Deta(ksi,eta,deltaPt5,deltaPt6) 0 funDN4Deta(ksi,eta,deltaPt6,deltaPt7) 0 ...
                        funDN5Deta(ksi,eta,deltaPt5) 0 funDN6Deta(ksi,eta,deltaPt6) 0 ...
                        funDN7Deta(ksi,eta,deltaPt7) 0 funDN8Deta(ksi,eta,deltaPt8) 0 ...
                        H1*funDN1Deta(ksi,eta,deltaPt7,deltaPt8) 0 H2*funDN2Deta(ksi,eta,deltaPt8,deltaPt5) 0 ...
                        H3*funDN3Deta(ksi,eta,deltaPt5,deltaPt6) 0 H4*funDN4Deta(ksi,eta,deltaPt6,deltaPt7) 0 ...
                        H5*funDN5Deta(ksi,eta,deltaPt5) 0 H6*funDN6Deta(ksi,eta,deltaPt6) 0 ...
                        H7*funDN7Deta(ksi,eta,deltaPt7) 0 H8*funDN8Deta(ksi,eta,deltaPt8) 0;
                        0 funDN1Dksi(ksi,eta,deltaPt7,deltaPt8) 0 funDN2Dksi(ksi,eta,deltaPt8,deltaPt5) ...
                        0 funDN3Dksi(ksi,eta,deltaPt5,deltaPt6) 0 funDN4Dksi(ksi,eta,deltaPt6,deltaPt7) ...
                        0 funDN5Dksi(ksi,eta,deltaPt5) 0 funDN6Dksi(ksi,eta,deltaPt6) ...
                        0 funDN7Dksi(ksi,eta,deltaPt7) 0 funDN8Dksi(ksi,eta,deltaPt8) ...
                        0 H1*funDN1Dksi(ksi,eta,deltaPt7,deltaPt8) 0 H2*funDN2Dksi(ksi,eta,deltaPt8,deltaPt5) ...
                        0 H3*funDN3Dksi(ksi,eta,deltaPt5,deltaPt6) 0 H4*funDN4Dksi(ksi,eta,deltaPt6,deltaPt7) ...
                        0 H5*funDN5Dksi(ksi,eta,deltaPt5) 0 H6*funDN6Dksi(ksi,eta,deltaPt6) ...
                        0 H7*funDN7Dksi(ksi,eta,deltaPt7) 0 H8*funDN8Dksi(ksi,eta,deltaPt8);
                        0 funDN1Deta(ksi,eta,deltaPt7,deltaPt8) 0 funDN2Deta(ksi,eta,deltaPt8,deltaPt5) ...
                        0 funDN3Deta(ksi,eta,deltaPt5,deltaPt6) 0 funDN4Deta(ksi,eta,deltaPt6,deltaPt7) ...
                        0 funDN5Deta(ksi,eta,deltaPt5) 0 funDN6Deta(ksi,eta,deltaPt6) ...
                        0 funDN7Deta(ksi,eta,deltaPt7) 0 funDN8Deta(ksi,eta,deltaPt8) ...
                        0 H1*funDN1Deta(ksi,eta,deltaPt7,deltaPt8) 0 H2*funDN2Deta(ksi,eta,deltaPt8,deltaPt5) ...
                        0 H3*funDN3Deta(ksi,eta,deltaPt5,deltaPt6) 0 H4*funDN4Deta(ksi,eta,deltaPt6,deltaPt7) ...
                        0 H5*funDN5Deta(ksi,eta,deltaPt5) 0 H6*funDN6Deta(ksi,eta,deltaPt6) ...
                        0 H7*funDN7Deta(ksi,eta,deltaPt7) 0 H8*funDN8Deta(ksi,eta,deltaPt8);];
                    
                elseif CrackOrNot == 1 && EleCrackHOrNot > 0 && EleCrackTipOrNot > 0  % Crack with crack tip
                    DN = [InvJ zeros(2,2); zeros(2,2) InvJ] * ...
                        [funDN1Dksi(ksi,eta,deltaPt7,deltaPt8) 0 funDN2Dksi(ksi,eta,deltaPt8,deltaPt5) 0 ...
                        funDN3Dksi(ksi,eta,deltaPt5,deltaPt6) 0 funDN4Dksi(ksi,eta,deltaPt6,deltaPt7) 0 ...
                        funDN5Dksi(ksi,eta,deltaPt5) 0 funDN6Dksi(ksi,eta,deltaPt6) 0 ...
                        funDN7Dksi(ksi,eta,deltaPt7) 0 funDN8Dksi(ksi,eta,deltaPt8) 0 ...
                        H1*funDN1Dksi(ksi,eta,deltaPt7,deltaPt8) 0 H2*funDN2Dksi(ksi,eta,deltaPt8,deltaPt5) 0 ...
                        H3*funDN3Dksi(ksi,eta,deltaPt5,deltaPt6) 0 H4*funDN4Dksi(ksi,eta,deltaPt6,deltaPt7) 0 ...
                        H5*funDN5Dksi(ksi,eta,deltaPt5) 0 H6*funDN6Dksi(ksi,eta,deltaPt6) 0 ...
                        H7*funDN7Dksi(ksi,eta,deltaPt7) 0 H8*funDN8Dksi(ksi,eta,deltaPt8) 0 ...
                        (C1*funDN1Dksi(ksi,eta,deltaPt7,deltaPt8)+N1*funDC1Dksi(ksi,eta,CrackTip))*deltaCNode1 0  (C2*funDN1Dksi(ksi,eta,deltaPt7,deltaPt8)+N1*funDC2Dksi(ksi,eta,CrackTip))*deltaCNode1 0 ...
                        (C3*funDN1Dksi(ksi,eta,deltaPt7,deltaPt8)+N1*funDC3Dksi(ksi,eta,CrackTip))*deltaCNode1 0  (C4*funDN1Dksi(ksi,eta,deltaPt7,deltaPt8)+N1*funDC4Dksi(ksi,eta,CrackTip))*deltaCNode1 0 ...
                        (C1*funDN2Dksi(ksi,eta,deltaPt8,deltaPt5)+N2*funDC1Dksi(ksi,eta,CrackTip))*deltaCNode2 0  (C2*funDN2Dksi(ksi,eta,deltaPt8,deltaPt5)+N2*funDC2Dksi(ksi,eta,CrackTip))*deltaCNode2 0 ...
                        (C3*funDN2Dksi(ksi,eta,deltaPt8,deltaPt5)+N2*funDC3Dksi(ksi,eta,CrackTip))*deltaCNode2 0  (C4*funDN2Dksi(ksi,eta,deltaPt8,deltaPt5)+N2*funDC4Dksi(ksi,eta,CrackTip))*deltaCNode2 0 ...
                        (C1*funDN3Dksi(ksi,eta,deltaPt5,deltaPt6)+N3*funDC1Dksi(ksi,eta,CrackTip))*deltaCNode3 0  (C2*funDN3Dksi(ksi,eta,deltaPt5,deltaPt6)+N3*funDC2Dksi(ksi,eta,CrackTip))*deltaCNode3 0 ...
                        (C3*funDN3Dksi(ksi,eta,deltaPt5,deltaPt6)+N3*funDC3Dksi(ksi,eta,CrackTip))*deltaCNode3 0  (C4*funDN3Dksi(ksi,eta,deltaPt5,deltaPt6)+N3*funDC4Dksi(ksi,eta,CrackTip))*deltaCNode3 0 ...
                        (C1*funDN4Dksi(ksi,eta,deltaPt6,deltaPt7)+N4*funDC1Dksi(ksi,eta,CrackTip))*deltaCNode4 0  (C2*funDN4Dksi(ksi,eta,deltaPt6,deltaPt7)+N4*funDC2Dksi(ksi,eta,CrackTip))*deltaCNode4 0 ...
                        (C3*funDN4Dksi(ksi,eta,deltaPt6,deltaPt7)+N4*funDC3Dksi(ksi,eta,CrackTip))*deltaCNode4 0  (C4*funDN4Dksi(ksi,eta,deltaPt6,deltaPt7)+N4*funDC4Dksi(ksi,eta,CrackTip))*deltaCNode4 0 ...
                        (C1*funDN5Dksi(ksi,eta,deltaPt5)+N5*funDC1Dksi(ksi,eta,CrackTip))*deltaCNode5 0 (C2*funDN5Dksi(ksi,eta,deltaPt5)+N5*funDC2Dksi(ksi,eta,CrackTip))*deltaCNode5 0 ...
                        (C3*funDN5Dksi(ksi,eta,deltaPt5)+N5*funDC3Dksi(ksi,eta,CrackTip))*deltaCNode5 0 (C4*funDN5Dksi(ksi,eta,deltaPt5)+N5*funDC4Dksi(ksi,eta,CrackTip))*deltaCNode5 0 ...
                        (C1*funDN6Dksi(ksi,eta,deltaPt6)+N6*funDC1Dksi(ksi,eta,CrackTip))*deltaCNode6 0 (C2*funDN6Dksi(ksi,eta,deltaPt6)+N6*funDC2Dksi(ksi,eta,CrackTip))*deltaCNode6 0 ...
                        (C3*funDN6Dksi(ksi,eta,deltaPt6)+N6*funDC3Dksi(ksi,eta,CrackTip))*deltaCNode6 0 (C4*funDN6Dksi(ksi,eta,deltaPt6)+N6*funDC4Dksi(ksi,eta,CrackTip))*deltaCNode6 0 ...
                        (C1*funDN7Dksi(ksi,eta,deltaPt7)+N7*funDC1Dksi(ksi,eta,CrackTip))*deltaCNode7 0 (C2*funDN7Dksi(ksi,eta,deltaPt7)+N7*funDC2Dksi(ksi,eta,CrackTip))*deltaCNode7 0 ...
                        (C3*funDN7Dksi(ksi,eta,deltaPt7)+N7*funDC3Dksi(ksi,eta,CrackTip))*deltaCNode7 0 (C4*funDN7Dksi(ksi,eta,deltaPt7)+N7*funDC4Dksi(ksi,eta,CrackTip))*deltaCNode7 0 ...
                        (C1*funDN8Dksi(ksi,eta,deltaPt8)+N8*funDC1Dksi(ksi,eta,CrackTip))*deltaCNode8 0 (C2*funDN8Dksi(ksi,eta,deltaPt8)+N8*funDC2Dksi(ksi,eta,CrackTip))*deltaCNode8 0 ...
                        (C3*funDN8Dksi(ksi,eta,deltaPt8)+N8*funDC3Dksi(ksi,eta,CrackTip))*deltaCNode8 0 (C4*funDN8Dksi(ksi,eta,deltaPt8)+N8*funDC4Dksi(ksi,eta,CrackTip))*deltaCNode8 0 ;
                        funDN1Deta(ksi,eta,deltaPt7,deltaPt8) 0 funDN2Deta(ksi,eta,deltaPt8,deltaPt5) 0 ...
                        funDN3Deta(ksi,eta,deltaPt5,deltaPt6) 0 funDN4Deta(ksi,eta,deltaPt6,deltaPt7) 0 ...
                        funDN5Deta(ksi,eta,deltaPt5) 0 funDN6Deta(ksi,eta,deltaPt6) 0 ...
                        funDN7Deta(ksi,eta,deltaPt7) 0 funDN8Deta(ksi,eta,deltaPt8) 0 ...
                        H1*funDN1Deta(ksi,eta,deltaPt7,deltaPt8) 0 H2*funDN2Deta(ksi,eta,deltaPt8,deltaPt5)   0 ...
                        H3*funDN3Deta(ksi,eta,deltaPt5,deltaPt6) 0 H4*funDN4Deta(ksi,eta,deltaPt6,deltaPt7)   0 ...
                        H5*funDN5Deta(ksi,eta,deltaPt5) 0 H6*funDN6Deta(ksi,eta,deltaPt6) 0 ...
                        H7*funDN7Deta(ksi,eta,deltaPt7) 0 H8*funDN8Deta(ksi,eta,deltaPt8) 0 ...
                        (C1*funDN1Deta(ksi,eta,deltaPt7,deltaPt8)+N1*funDC1Deta(ksi,eta,CrackTip))*deltaCNode1 0  (C2*funDN1Deta(ksi,eta,deltaPt7,deltaPt8)+N1*funDC2Deta(ksi,eta,CrackTip))*deltaCNode1 0 ...
                        (C3*funDN1Deta(ksi,eta,deltaPt7,deltaPt8)+N1*funDC3Deta(ksi,eta,CrackTip))*deltaCNode1 0  (C4*funDN1Deta(ksi,eta,deltaPt7,deltaPt8)+N1*funDC4Deta(ksi,eta,CrackTip))*deltaCNode1 0 ...
                        (C1*funDN2Deta(ksi,eta,deltaPt8,deltaPt5)+N2*funDC1Deta(ksi,eta,CrackTip))*deltaCNode2 0  (C2*funDN2Deta(ksi,eta,deltaPt8,deltaPt5)+N2*funDC2Deta(ksi,eta,CrackTip))*deltaCNode2 0 ...
                        (C3*funDN2Deta(ksi,eta,deltaPt8,deltaPt5)+N2*funDC3Deta(ksi,eta,CrackTip))*deltaCNode2 0  (C4*funDN2Deta(ksi,eta,deltaPt8,deltaPt5)+N2*funDC4Deta(ksi,eta,CrackTip))*deltaCNode2 0 ...
                        (C1*funDN3Deta(ksi,eta,deltaPt5,deltaPt6)+N3*funDC1Deta(ksi,eta,CrackTip))*deltaCNode3 0  (C2*funDN3Deta(ksi,eta,deltaPt5,deltaPt6)+N3*funDC2Deta(ksi,eta,CrackTip))*deltaCNode3 0 ...
                        (C3*funDN3Deta(ksi,eta,deltaPt5,deltaPt6)+N3*funDC3Deta(ksi,eta,CrackTip))*deltaCNode3 0  (C4*funDN3Deta(ksi,eta,deltaPt5,deltaPt6)+N3*funDC4Deta(ksi,eta,CrackTip))*deltaCNode3 0 ...
                        (C1*funDN4Deta(ksi,eta,deltaPt6,deltaPt7)+N4*funDC1Deta(ksi,eta,CrackTip))*deltaCNode4 0  (C2*funDN4Deta(ksi,eta,deltaPt6,deltaPt7)+N4*funDC2Deta(ksi,eta,CrackTip))*deltaCNode4 0 ...
                        (C3*funDN4Deta(ksi,eta,deltaPt6,deltaPt7)+N4*funDC3Deta(ksi,eta,CrackTip))*deltaCNode4 0  (C4*funDN4Deta(ksi,eta,deltaPt6,deltaPt7)+N4*funDC4Deta(ksi,eta,CrackTip))*deltaCNode4 0 ...
                        (C1*funDN5Deta(ksi,eta,deltaPt5)+N5*funDC1Deta(ksi,eta,CrackTip))*deltaCNode5 0 (C2*funDN5Deta(ksi,eta,deltaPt5)+N5*funDC2Deta(ksi,eta,CrackTip))*deltaCNode5 0 ...
                        (C3*funDN5Deta(ksi,eta,deltaPt5)+N5*funDC3Deta(ksi,eta,CrackTip))*deltaCNode5 0 (C4*funDN5Deta(ksi,eta,deltaPt5)+N5*funDC4Deta(ksi,eta,CrackTip))*deltaCNode5 0 ...
                        (C1*funDN6Deta(ksi,eta,deltaPt6)+N6*funDC1Deta(ksi,eta,CrackTip))*deltaCNode6 0 (C2*funDN6Deta(ksi,eta,deltaPt6)+N6*funDC2Deta(ksi,eta,CrackTip))*deltaCNode6 0 ...
                        (C3*funDN6Deta(ksi,eta,deltaPt6)+N6*funDC3Deta(ksi,eta,CrackTip))*deltaCNode6 0 (C4*funDN6Deta(ksi,eta,deltaPt6)+N6*funDC4Deta(ksi,eta,CrackTip))*deltaCNode6 0 ...
                        (C1*funDN7Deta(ksi,eta,deltaPt7)+N7*funDC1Deta(ksi,eta,CrackTip))*deltaCNode7 0 (C2*funDN7Deta(ksi,eta,deltaPt7)+N7*funDC2Deta(ksi,eta,CrackTip))*deltaCNode7 0 ...
                        (C3*funDN7Deta(ksi,eta,deltaPt7)+N7*funDC3Deta(ksi,eta,CrackTip))*deltaCNode7 0 (C4*funDN7Deta(ksi,eta,deltaPt7)+N7*funDC4Deta(ksi,eta,CrackTip))*deltaCNode7 0 ...
                        (C1*funDN8Deta(ksi,eta,deltaPt8)+N8*funDC1Deta(ksi,eta,CrackTip))*deltaCNode8 0 (C2*funDN8Deta(ksi,eta,deltaPt8)+N8*funDC2Deta(ksi,eta,CrackTip))*deltaCNode8 0 ...
                        (C3*funDN8Deta(ksi,eta,deltaPt8)+N8*funDC3Deta(ksi,eta,CrackTip))*deltaCNode8 0 (C4*funDN8Deta(ksi,eta,deltaPt8)+N8*funDC4Deta(ksi,eta,CrackTip))*deltaCNode8 0 ;
                        0 funDN1Dksi(ksi,eta,deltaPt7,deltaPt8) 0 funDN2Dksi(ksi,eta,deltaPt8,deltaPt5) ...
                        0 funDN3Dksi(ksi,eta,deltaPt5,deltaPt6) 0 funDN4Dksi(ksi,eta,deltaPt6,deltaPt7) ...
                        0 funDN5Dksi(ksi,eta,deltaPt5) 0 funDN6Dksi(ksi,eta,deltaPt6) ...
                        0 funDN7Dksi(ksi,eta,deltaPt7) 0 funDN8Dksi(ksi,eta,deltaPt8) ...
                        0 H1*funDN1Dksi(ksi,eta,deltaPt7,deltaPt8) 0 H2*funDN2Dksi(ksi,eta,deltaPt8,deltaPt5) ...
                        0 H3*funDN3Dksi(ksi,eta,deltaPt5,deltaPt6) 0 H4*funDN4Dksi(ksi,eta,deltaPt6,deltaPt7) ...
                        0 H5*funDN5Dksi(ksi,eta,deltaPt5) 0 H6*funDN6Dksi(ksi,eta,deltaPt6) ...
                        0 H7*funDN7Dksi(ksi,eta,deltaPt7) 0 H8*funDN8Dksi(ksi,eta,deltaPt8) ...
                        0 (C1*funDN1Dksi(ksi,eta,deltaPt7,deltaPt8)+N1*funDC1Dksi(ksi,eta,CrackTip))*deltaCNode1 0  (C2*funDN1Dksi(ksi,eta,deltaPt7,deltaPt8)+N1*funDC2Dksi(ksi,eta,CrackTip))*deltaCNode1 ...
                        0 (C3*funDN1Dksi(ksi,eta,deltaPt7,deltaPt8)+N1*funDC3Dksi(ksi,eta,CrackTip))*deltaCNode1 0  (C4*funDN1Dksi(ksi,eta,deltaPt7,deltaPt8)+N1*funDC4Dksi(ksi,eta,CrackTip))*deltaCNode1 ...
                        0 (C1*funDN2Dksi(ksi,eta,deltaPt8,deltaPt5)+N2*funDC1Dksi(ksi,eta,CrackTip))*deltaCNode2 0  (C2*funDN2Dksi(ksi,eta,deltaPt8,deltaPt5)+N2*funDC2Dksi(ksi,eta,CrackTip))*deltaCNode2 ...
                        0 (C3*funDN2Dksi(ksi,eta,deltaPt8,deltaPt5)+N2*funDC3Dksi(ksi,eta,CrackTip))*deltaCNode2 0  (C4*funDN2Dksi(ksi,eta,deltaPt8,deltaPt5)+N2*funDC4Dksi(ksi,eta,CrackTip))*deltaCNode2 ...
                        0 (C1*funDN3Dksi(ksi,eta,deltaPt5,deltaPt6)+N3*funDC1Dksi(ksi,eta,CrackTip))*deltaCNode3 0  (C2*funDN3Dksi(ksi,eta,deltaPt5,deltaPt6)+N3*funDC2Dksi(ksi,eta,CrackTip))*deltaCNode3 ...
                        0 (C3*funDN3Dksi(ksi,eta,deltaPt5,deltaPt6)+N3*funDC3Dksi(ksi,eta,CrackTip))*deltaCNode3 0  (C4*funDN3Dksi(ksi,eta,deltaPt5,deltaPt6)+N3*funDC4Dksi(ksi,eta,CrackTip))*deltaCNode3 ...
                        0 (C1*funDN4Dksi(ksi,eta,deltaPt6,deltaPt7)+N4*funDC1Dksi(ksi,eta,CrackTip))*deltaCNode4 0  (C2*funDN4Dksi(ksi,eta,deltaPt6,deltaPt7)+N4*funDC2Dksi(ksi,eta,CrackTip))*deltaCNode4 ...
                        0 (C3*funDN4Dksi(ksi,eta,deltaPt6,deltaPt7)+N4*funDC3Dksi(ksi,eta,CrackTip))*deltaCNode4 0  (C4*funDN4Dksi(ksi,eta,deltaPt6,deltaPt7)+N4*funDC4Dksi(ksi,eta,CrackTip))*deltaCNode4 ...
                        0 (C1*funDN5Dksi(ksi,eta,deltaPt5)+N5*funDC1Dksi(ksi,eta,CrackTip))*deltaCNode5 0 (C2*funDN5Dksi(ksi,eta,deltaPt5)+N5*funDC2Dksi(ksi,eta,CrackTip))*deltaCNode5 ...
                        0 (C3*funDN5Dksi(ksi,eta,deltaPt5)+N5*funDC3Dksi(ksi,eta,CrackTip))*deltaCNode5 0 (C4*funDN5Dksi(ksi,eta,deltaPt5)+N5*funDC4Dksi(ksi,eta,CrackTip))*deltaCNode5 ...
                        0 (C1*funDN6Dksi(ksi,eta,deltaPt6)+N6*funDC1Dksi(ksi,eta,CrackTip))*deltaCNode6 0 (C2*funDN6Dksi(ksi,eta,deltaPt6)+N6*funDC2Dksi(ksi,eta,CrackTip))*deltaCNode6 ...
                        0 (C3*funDN6Dksi(ksi,eta,deltaPt6)+N6*funDC3Dksi(ksi,eta,CrackTip))*deltaCNode6 0 (C4*funDN6Dksi(ksi,eta,deltaPt6)+N6*funDC4Dksi(ksi,eta,CrackTip))*deltaCNode6 ...
                        0 (C1*funDN7Dksi(ksi,eta,deltaPt7)+N7*funDC1Dksi(ksi,eta,CrackTip))*deltaCNode7 0 (C2*funDN7Dksi(ksi,eta,deltaPt7)+N7*funDC2Dksi(ksi,eta,CrackTip))*deltaCNode7 ...
                        0 (C3*funDN7Dksi(ksi,eta,deltaPt7)+N7*funDC3Dksi(ksi,eta,CrackTip))*deltaCNode7 0 (C4*funDN7Dksi(ksi,eta,deltaPt7)+N7*funDC4Dksi(ksi,eta,CrackTip))*deltaCNode7 ...
                        0 (C1*funDN8Dksi(ksi,eta,deltaPt8)+N8*funDC1Dksi(ksi,eta,CrackTip))*deltaCNode8 0 (C2*funDN8Dksi(ksi,eta,deltaPt8)+N8*funDC2Dksi(ksi,eta,CrackTip))*deltaCNode8 ...
                        0 (C3*funDN8Dksi(ksi,eta,deltaPt8)+N8*funDC3Dksi(ksi,eta,CrackTip))*deltaCNode8 0 (C4*funDN8Dksi(ksi,eta,deltaPt8)+N8*funDC4Dksi(ksi,eta,CrackTip))*deltaCNode8 ;
                        0 funDN1Deta(ksi,eta,deltaPt7,deltaPt8) 0 funDN2Deta(ksi,eta,deltaPt8,deltaPt5) ...
                        0 funDN3Deta(ksi,eta,deltaPt5,deltaPt6) 0 funDN4Deta(ksi,eta,deltaPt6,deltaPt7) ...
                        0 funDN5Deta(ksi,eta,deltaPt5) 0 funDN6Deta(ksi,eta,deltaPt6) ...
                        0 funDN7Deta(ksi,eta,deltaPt7) 0 funDN8Deta(ksi,eta,deltaPt8) ...
                        0 H1*funDN1Deta(ksi,eta,deltaPt7,deltaPt8) 0 H2*funDN2Deta(ksi,eta,deltaPt8,deltaPt5) ...
                        0 H3*funDN3Deta(ksi,eta,deltaPt5,deltaPt6) 0 H4*funDN4Deta(ksi,eta,deltaPt6,deltaPt7) ...
                        0 H5*funDN5Deta(ksi,eta,deltaPt5) 0 H6*funDN6Deta(ksi,eta,deltaPt6) ...
                        0 H7*funDN7Deta(ksi,eta,deltaPt7) 0 H8*funDN8Deta(ksi,eta,deltaPt8) ...
                        0 (C1*funDN1Deta(ksi,eta,deltaPt7,deltaPt8)+N1*funDC1Deta(ksi,eta,CrackTip))*deltaCNode1 0  (C2*funDN1Deta(ksi,eta,deltaPt7,deltaPt8)+N1*funDC2Deta(ksi,eta,CrackTip))*deltaCNode1 ...
                        0 (C3*funDN1Deta(ksi,eta,deltaPt7,deltaPt8)+N1*funDC3Deta(ksi,eta,CrackTip))*deltaCNode1 0  (C4*funDN1Deta(ksi,eta,deltaPt7,deltaPt8)+N1*funDC4Deta(ksi,eta,CrackTip))*deltaCNode1 ...
                        0 (C1*funDN2Deta(ksi,eta,deltaPt8,deltaPt5)+N2*funDC1Deta(ksi,eta,CrackTip))*deltaCNode2 0  (C2*funDN2Deta(ksi,eta,deltaPt8,deltaPt5)+N2*funDC2Deta(ksi,eta,CrackTip))*deltaCNode2 ...
                        0 (C3*funDN2Deta(ksi,eta,deltaPt8,deltaPt5)+N2*funDC3Deta(ksi,eta,CrackTip))*deltaCNode2 0  (C4*funDN2Deta(ksi,eta,deltaPt8,deltaPt5)+N2*funDC4Deta(ksi,eta,CrackTip))*deltaCNode2 ...
                        0 (C1*funDN3Deta(ksi,eta,deltaPt5,deltaPt6)+N3*funDC1Deta(ksi,eta,CrackTip))*deltaCNode3 0  (C2*funDN3Deta(ksi,eta,deltaPt5,deltaPt6)+N3*funDC2Deta(ksi,eta,CrackTip))*deltaCNode3 ...
                        0 (C3*funDN3Deta(ksi,eta,deltaPt5,deltaPt6)+N3*funDC3Deta(ksi,eta,CrackTip))*deltaCNode3 0  (C4*funDN3Deta(ksi,eta,deltaPt5,deltaPt6)+N3*funDC4Deta(ksi,eta,CrackTip))*deltaCNode3 ...
                        0 (C1*funDN4Deta(ksi,eta,deltaPt6,deltaPt7)+N4*funDC1Deta(ksi,eta,CrackTip))*deltaCNode4 0  (C2*funDN4Deta(ksi,eta,deltaPt6,deltaPt7)+N4*funDC2Deta(ksi,eta,CrackTip))*deltaCNode4 ...
                        0 (C3*funDN4Deta(ksi,eta,deltaPt6,deltaPt7)+N4*funDC3Deta(ksi,eta,CrackTip))*deltaCNode4 0  (C4*funDN4Deta(ksi,eta,deltaPt6,deltaPt7)+N4*funDC4Deta(ksi,eta,CrackTip))*deltaCNode4 ...
                        0 (C1*funDN5Deta(ksi,eta,deltaPt5)+N5*funDC1Deta(ksi,eta,CrackTip))*deltaCNode5 0 (C2*funDN5Deta(ksi,eta,deltaPt5)+N5*funDC2Deta(ksi,eta,CrackTip))*deltaCNode5 ...
                        0 (C3*funDN5Deta(ksi,eta,deltaPt5)+N5*funDC3Deta(ksi,eta,CrackTip))*deltaCNode5 0 (C4*funDN5Deta(ksi,eta,deltaPt5)+N5*funDC4Deta(ksi,eta,CrackTip))*deltaCNode5 ...
                        0 (C1*funDN6Deta(ksi,eta,deltaPt6)+N6*funDC1Deta(ksi,eta,CrackTip))*deltaCNode6 0 (C2*funDN6Deta(ksi,eta,deltaPt6)+N6*funDC2Deta(ksi,eta,CrackTip))*deltaCNode6 ...
                        0 (C3*funDN6Deta(ksi,eta,deltaPt6)+N6*funDC3Deta(ksi,eta,CrackTip))*deltaCNode6 0 (C4*funDN6Deta(ksi,eta,deltaPt6)+N6*funDC4Deta(ksi,eta,CrackTip))*deltaCNode6 ...
                        0 (C1*funDN7Deta(ksi,eta,deltaPt7)+N7*funDC1Deta(ksi,eta,CrackTip))*deltaCNode7 0 (C2*funDN7Deta(ksi,eta,deltaPt7)+N7*funDC2Deta(ksi,eta,CrackTip))*deltaCNode7 ...
                        0 (C3*funDN7Deta(ksi,eta,deltaPt7)+N7*funDC3Deta(ksi,eta,CrackTip))*deltaCNode7 0 (C4*funDN7Deta(ksi,eta,deltaPt7)+N7*funDC4Deta(ksi,eta,CrackTip))*deltaCNode7 ...
                        0 (C1*funDN8Deta(ksi,eta,deltaPt8)+N8*funDC1Deta(ksi,eta,CrackTip))*deltaCNode8 0 (C2*funDN8Deta(ksi,eta,deltaPt8)+N8*funDC2Deta(ksi,eta,CrackTip))*deltaCNode8 ...
                        0 (C3*funDN8Deta(ksi,eta,deltaPt8)+N8*funDC3Deta(ksi,eta,CrackTip))*deltaCNode8 0 (C4*funDN8Deta(ksi,eta,deltaPt8)+N8*funDC4Deta(ksi,eta,CrackTip))*deltaCNode8 ; ];
                          
                end
                
                % ------- Comment: Calculate DivFMinusW ---------
                % tempDUDX = DN*FMinusW1(temp);
                % DivFMinusW1 = tempDUDX(1)+tempDUDX(4);
                % tempDUDX = DN*FMinusW2(temp);
                % DivFMinusW2 = tempDUDX(1)+tempDUDX(4);
                % DivFMinusW = [DivFMinusW1;DivFMinusW2;DivFMinusW1;DivFMinusW2;
                %               DivFMinusW1;DivFMinusW2;DivFMinusW1;DivFMinusW2;
                %               DivFMinusW1;DivFMinusW2;DivFMinusW1;DivFMinusW2;
                %               DivFMinusW1;DivFMinusW2;DivFMinusW1;DivFMinusW2];
                % ------ End of comment of calculating DivFMinusW ------
                
                % ------ Construct A matrix ------
                tempA = tempA + Jacobian*weightksi*weighteta * ( (beta+alpha)*(DN')*(DN) + mu*(N')*N  );
                
                % ------ Construct b vector ------
                tempb = tempb + Jacobian*weightksi*weighteta * ( beta*diag((DN')*(FMinusW(tempIndexF)')) + mu*N'*N*(UMinusv(tempIndexU)) +(alpha)*(DN')*DN*U(tempIndexU) );
                  
            %end
        %end
        end
        
        
        % Store tempA and tempb
        if IterStep==1
            [IndexAXX,IndexAYY] = ndgrid(tempIndexU,tempIndexU);
            INDEXAIpar{j}=IndexAXX(:); INDEXAJpar{j}=IndexAYY(:); INDEXAVALpar{j}=tempA(:);
        end
        INDEXBIpar{j}=tempIndexU(:); INDEXBVALpar{j}=tempb(:);   
        
    end
    
    close(hbar);
    
    
    % ====== Initialize A matrix and b vector ======
    A = sparse(FEMSize+DIM*NodesNumPerEle, FEMSize+DIM*NodesNumPerEle);   b=sparse(FEMSize+DIM*NodesNumPerEle,1);
    for eleInd = 1:size(elementsFEM,1)
        A = A + sparse(INDEXAIpar{eleInd}, INDEXAJpar{eleInd}, INDEXAVALpar{eleInd}, FEMSize+DIM*NodesNumPerEle , FEMSize+DIM*NodesNumPerEle) ;
        b = b + sparse(INDEXBIpar{eleInd},ones(length(INDEXBIpar{eleInd}),1),INDEXBVALpar{eleInd}, FEMSize+DIM*NodesNumPerEle , 1) ;
    end
% "DIM*NodesNumPerEle" is because I put all the zeros in elementsFEM to the end, and in NC function there are "DIM*NodesNumPerEle" entries




    %if IterStep == 1
%        AMatrixRegularized = A + 1e-3*max(diag(A))*speye(size(A,1),size(A,2));
        % AMatrixRegularized = A; tempVal = 1e-3*max(diag(A));
        % for tempi = 1:size(A,1)
        %     AMatrixRegularized(tempi,tempi) = AMatrixRegularized(tempi,tempi) + tempVal;
        % send
    %end
    
    coordsIndexInvolved = unique(elementsFEM(:,1:8)); % Need modification for triangle elementsFEM
    if CrackOrNot == 1
        coordsIndexInvolvedH  = EnrHAndTipEleIndex(:,2);
        if CrackTipOrNot == 0
            coordsIndexInvolvedC1 = EnrTipEleIndex(:,2); coordsIndexInvolvedC2 = EnrTipEleIndex(:,2)+1;
            coordsIndexInvolvedC3 = EnrTipEleIndex(:,2)+2; coordsIndexInvolvedC4 = EnrTipEleIndex(:,2)+3;
            coordsIndexInvolved   = [coordsIndexInvolved;coordsIndexInvolvedH;coordsIndexInvolvedC1;coordsIndexInvolvedC2;coordsIndexInvolvedC3;coordsIndexInvolvedC4];
        else
            coordsIndexInvolved   = [coordsIndexInvolved;coordsIndexInvolvedH];
        end
    end
    UIndexInvolved = zeros(2*(length(coordsIndexInvolved)-1),1);
    % Not including the first 0-th entry
    for tempi = 1:(size(coordsIndexInvolved,1)-1)
        UIndexInvolved(2*tempi-1:2*tempi) = [2*coordsIndexInvolved(tempi+1)-1; 2*coordsIndexInvolved(tempi+1)];
    end
    
    % ========= Set Dirichlet and Neumann boundary conditions =========
    if isempty(dirichlet) ~= 1
        dirichlettemp = [2*dirichlet(:); 2*dirichlet(:)-1];
    else
        dirichlettemp = [];
    end
%     if isempty(neumann) ~= 1
%         neumanntemp = [2*neumann(:,1); 2*neumann(:,1)-1; 2*neumann(:,2); 2*neumann(:,2)-1];
%     else
%         neumanntemp = [];
%     end
    FreeNodes = setdiff(UIndexInvolved,unique([dirichlettemp]));
    
    % ========= Neumann conditions ===========
    % Last step boundary condition force
    BCForce = - 1/winstepsize * FSubpb1;
    for tempj = 1:size(neumann,1)
        b(2*neumann(tempj,1:2)-1) = b(2*neumann(tempj,1:2)-1) + 0.5*norm(coordinatesFEM(neumann(tempj,1),:)-coordinatesFEM(neumann(tempj,2),:)) ...
          *( ( BCForce(4*neumann(tempj,1:2)-3) * neumann(tempj,3) + BCForce(4*neumann(tempj,1:2)-1) * neumann(tempj,4) ) );
        b(2*neumann(tempj,1:2))   = b(2*neumann(tempj,1:2)) + 0.5*norm(coordinatesFEM(neumann(tempj,1),:)-coordinatesFEM(neumann(tempj,2),:)) ...
         * (( BCForce(4*neumann(tempj,1:2)-2) * neumann(tempj,3) + BCForce(4*neumann(tempj,1:2)) * neumann(tempj,4) ) );
    end
    
    % ========= Dirichlet conditions ==========
    UhatOld = Uhat; Uhat = sparse(FEMSize+DIM*NodesNumPerEle,1);
    for tempi = 1:DIM
        Uhat(DIM*unique(dirichlet)-(tempi-1)) = U(DIM*unique(dirichlet)-(tempi-1)); 
    end
    b = b - A * Uhat;
     
    % ========= Solve FEM problem ===========
    Uhat(FreeNodes) = A(FreeNodes,FreeNodes) \ b(FreeNodes);
    UhatNew = Uhat(1:end-DIM*NodesNumPerEle);
     
    if norm(UhatNew-UhatOld)/sqrt(length(UhatOld)) < 1e-3
        ConvergeOrNot = 1;
    end
    
    Uhat = UhatNew;
    
    
    
end



%% ========= subroutines for  FEM Q4 shape function derivatives ========

function DN5Dksi=funDN5Dksi(ksi,eta,deltaPt5)
%DN5Dksi = deltaPt5*0.5*(1-abs(eta)) ;
DN5Dksi = deltaPt5*0.5*(1-eta^2);

function DN5Deta=funDN5Deta(ksi,eta,deltaPt5)
%DN5Deta = deltaPt5*0.5*(1+ksi)*sign(-eta);
DN5Deta = deltaPt5*(-1)*(1+ksi)*eta;

function DN6Dksi=funDN6Dksi(ksi,eta,deltaPt6)
%DN6Dksi = deltaPt6*0.5*(1+eta)*sign(-ksi);
DN6Dksi = deltaPt6*(-1)*(1+eta)*ksi;

function DN6Deta=funDN6Deta(ksi,eta,deltaPt6)
% DN6Deta = deltaPt6*0.5*(1-abs(ksi));
DN6Deta = deltaPt6*0.5*(1-ksi^2);

function DN7Dksi=funDN7Dksi(ksi,eta,deltaPt7)
%DN7Dksi = deltaPt7*0.5*(-1)*(1-abs(eta));
DN7Dksi = deltaPt7*(-0.5)*(1-eta^2);

function DN7Deta=funDN7Deta(ksi,eta,deltaPt7)
% DN7Deta = deltaPt7*0.5*(1-ksi)*sign(-eta);
DN7Deta = deltaPt7*(-1)*(1-ksi)*eta;

function DN8Dksi=funDN8Dksi(ksi,eta,deltaPt8)
% DN8Dksi = deltaPt8*0.5*(1-eta)*sign(-ksi);
DN8Dksi = deltaPt8*(-1)*(1-eta)*ksi;

function DN8Deta=funDN8Deta(ksi,eta,deltaPt8)
% DN8Deta = deltaPt8*0.5*(-1)*(1-abs(ksi));
DN8Deta = deltaPt8*(-0.5)*(1-ksi^2);

function DN1Dksi = funDN1Dksi(ksi,eta,deltaPt7,deltaPt8)
DN1Dksi = -0.25*(1-eta)-0.5*( funDN7Dksi(ksi,eta,deltaPt7) + funDN8Dksi(ksi,eta,deltaPt8) );

function DN1Deta = funDN1Deta(ksi,eta,deltaPt7,deltaPt8)
DN1Deta = -0.25*(1-ksi)-0.5*( funDN7Deta(ksi,eta,deltaPt7) + funDN8Deta(ksi,eta,deltaPt8) );

function DN2Dksi = funDN2Dksi(ksi,eta,deltaPt8,deltaPt5)
DN2Dksi = 0.25*(1-eta)-0.5*( funDN8Dksi(ksi,eta,deltaPt8) + funDN5Dksi(ksi,eta,deltaPt5) );

function DN2Deta = funDN2Deta(ksi,eta,deltaPt8,deltaPt5)
DN2Deta = -0.25*(1+ksi)-0.5*( funDN8Deta(ksi,eta,deltaPt8) + funDN5Deta(ksi,eta,deltaPt5) );

function DN3Dksi = funDN3Dksi(ksi,eta,deltaPt5,deltaPt6)
DN3Dksi = 0.25*(1+eta)-0.5*( funDN5Dksi(ksi,eta,deltaPt5) + funDN6Dksi(ksi,eta,deltaPt6) );

function DN3Deta = funDN3Deta(ksi,eta,deltaPt5,deltaPt6)
DN3Deta = 0.25*(1+ksi)-0.5*( funDN5Deta(ksi,eta,deltaPt5) + funDN6Deta(ksi,eta,deltaPt6) );

function DN4Dksi = funDN4Dksi(ksi,eta,deltaPt6,deltaPt7)
DN4Dksi = -0.25*(1+eta)-0.5*( funDN6Dksi(ksi,eta,deltaPt6) + funDN7Dksi(ksi,eta,deltaPt7) );

function DN4Deta = funDN4Deta(ksi,eta,deltaPt6,deltaPt7)
DN4Deta = 0.25*(1-ksi)-0.5*( funDN6Deta(ksi,eta,deltaPt6) + funDN7Deta(ksi,eta,deltaPt7) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function DC1Dksi = funDC1Dksi(ksi,eta,CrackTip)
if (ksi^2+eta^2) > abs(eps)
    DC1Dksi = sqrt(ksi^2+eta^2) * cos(0.5*atan2((ksi-CrackTip(1)),(eta-CrackTip(2))))*0.5*(1/(1+((ksi-CrackTip(1))/(eta-CrackTip(2)))^2))*(1/(eta-CrackTip(2))) + ...
        (1/sqrt(ksi^2+eta^2))*ksi*sin(0.5*atan2((ksi-CrackTip(1)),(eta-CrackTip(2))));
else
    DC1Dksi = 0;
end

function DC1Deta = funDC1Deta(ksi,eta,CrackTip)
if (ksi^2+eta^2) > abs(eps)
    DC1Deta = sqrt(ksi^2+eta^2) * cos(0.5*atan2((ksi-CrackTip(1)),(eta-CrackTip(2))))*0.5*(1/(1+((ksi-CrackTip(1))/(eta-CrackTip(2)))^2))*(-(ksi-CrackTip(1))/(eta-CrackTip(2))^2) + ...
        (1/sqrt(ksi^2+eta^2))*eta*sin(0.5*atan2((ksi-CrackTip(1)),(eta-CrackTip(2))));
else
    DC1Deta = 0;
end

function DC2Dksi = funDC2Dksi(ksi,eta,CrackTip)
if (ksi^2+eta^2) > abs(eps)
    DC2Dksi = sqrt(ksi^2+eta^2) * (-sin(0.5*atan2((ksi-CrackTip(1)),(eta-CrackTip(2)))))*0.5*(1/(1+((ksi-CrackTip(1))/(eta-CrackTip(2)))^2))*(1/(eta-CrackTip(2))) + ...
        (1/sqrt(ksi^2+eta^2))*ksi*cos(0.5*atan2((ksi-CrackTip(1)),(eta-CrackTip(2))));
else
    DC2Dksi = 0;
end

function DC2Deta = funDC2Deta(ksi,eta,CrackTip)
if (ksi^2+eta^2) > abs(eps)
    DC2Deta = sqrt(ksi^2+eta^2) * -sin(0.5*atan2((ksi-CrackTip(1)),(eta-CrackTip(2))))*0.5*(1/(1+((ksi-CrackTip(1))/(eta-CrackTip(2)))^2))*(-(ksi-CrackTip(1))/(eta-CrackTip(2))^2) + ...
        (1/sqrt(ksi^2+eta^2))*eta*cos(0.5*atan2((ksi-CrackTip(1)),(eta-CrackTip(2))));
else
    DC2Deta = 0;
end
    
function DC3Dksi = funDC3Dksi(ksi,eta,CrackTip)
if (ksi^2+eta^2) > abs(eps)
    C1 = sqrt(ksi^2+eta^2) * sin(0.5*atan2((ksi-CrackTip(1)),(eta-CrackTip(2))));
    DC3Dksi = funDC1Dksi(ksi,eta,CrackTip)*sin(atan2((ksi-CrackTip(1)),(eta-CrackTip(2)))) + C1*cos(atan2((ksi-CrackTip(1)),(eta-CrackTip(2))))*(1/(1+((ksi-CrackTip(1))/(eta-CrackTip(2)))^2))*(1/(eta-CrackTip(2)));
else
    DC3Dksi = 0;
end

function DC3Deta = funDC3Deta(ksi,eta,CrackTip)
if (ksi^2+eta^2) > abs(eps)
    C1 = sqrt(ksi^2+eta^2) * sin(0.5*atan2((ksi-CrackTip(1)),(eta-CrackTip(2))));
    DC3Deta = funDC1Deta(ksi,eta,CrackTip)*sin(atan2((ksi-CrackTip(1)),(eta-CrackTip(2)))) + C1*cos(atan2((ksi-CrackTip(1)),(eta-CrackTip(2))))*(1/(1+((ksi-CrackTip(1))/(eta-CrackTip(2)))^2))*(-(ksi-CrackTip(1))/(eta-CrackTip(2))^2);
else
    DC3Deta = 0;
end
    
function DC4Dksi = funDC4Dksi(ksi,eta,CrackTip)
if (ksi^2+eta^2) > abs(eps)
    C2 = sqrt(ksi^2+eta^2) * cos(0.5*atan2((ksi-CrackTip(1)),(eta-CrackTip(2))));
    DC4Dksi = funDC2Dksi(ksi,eta,CrackTip)*sin(atan2((ksi-CrackTip(1)),(eta-CrackTip(2)))) + C2*(cos(atan2((ksi-CrackTip(1)),(eta-CrackTip(2)))))*(1/(1+((ksi-CrackTip(1))/(eta-CrackTip(2)))^2))*(1/(eta-CrackTip(2)));
else
    DC4Dksi = 0;
end

function DC4Deta = funDC4Deta(ksi,eta,CrackTip)
if (ksi^2+eta^2) > abs(eps)
    C2 = sqrt(ksi^2+eta^2) * cos(0.5*atan2((ksi-CrackTip(1)),(eta-CrackTip(2))));
    DC4Deta = funDC2Deta(ksi,eta,CrackTip)*sin(atan2((ksi-CrackTip(1)),(eta-CrackTip(2)))) + C2*(cos(atan2((ksi-CrackTip(1)),(eta-CrackTip(2)))))*(1/(1+((ksi-CrackTip(1))/(eta-CrackTip(2)))^2))*(-(ksi-CrackTip(1))/(eta-CrackTip(2))^2);
else
    DC4Deta = 0;
end




