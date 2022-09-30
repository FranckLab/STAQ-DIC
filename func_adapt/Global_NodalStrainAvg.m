function [Strain,StrainGaussPt,CoordGaussPt] = Global_NodalStrainAvg(coordinates,elements,U,GaussPtOrder, ...
    winstepsize,LevelNo,...
    EnrHAndTipEleIndex,EnrTipEleIndex,CrackOrNot,CrackPath1,CrackPath2,CrackTip,CrackTipOrNot)

DIM = 2; NodesNumPerEle = 4;

% ====== Info of CrackTip path-line ======
k_tip = -0.5*(CrackPath1(1)/CrackPath1(2)+CrackPath2(1)/CrackPath2(2));
CrackPathCen = [k_tip/(CrackTip(2)-k_tip*CrackTip(1)), -1/(CrackTip(2)-k_tip*CrackTip(1))];

CrackTipOrNot = 1-CrackTipOrNot; % Keep consistent with old version codes.
if CrackOrNot == 0
    FEMSize = size(coordinates,1);
elseif CrackOrNot>0 && CrackTipOrNot == 0
    FEMSize = size(coordinates,1) + size(EnrHAndTipEleIndex,1) + 4*size(EnrTipEleIndex,1);
else
    FEMSize = size(coordinates,1) + size(EnrHAndTipEleIndex,1);
end

% ====== Initialize strain Gauss points ======
Strain = 0; U = [U;zeros(DIM*NodesNumPerEle,1)];
FStrainAvgTimes = zeros(4*size(coordinates,1),1); FStrain = zeros(4*size(coordinates,1),1);
StrainGaussPt = zeros(4*size(elements,1),4); CoordGaussPt = zeros(4*size(elements,1),2);

% ====== Gaussian quadrature parameter ======
switch GaussPtOrder
    case 2 % ------ 2*2 Gauss points ------
        gqpt1 = -1/sqrt(3); gqpt2 = 1/sqrt(3); gqpt = [gqpt1,gqpt2]; 
        gqwt1 = 1; gqwt2 = 1; gqwt = [gqwt1,gqwt2];
    case 3 % ------ 3*3 Gauss points ------
        gqpt1 = 0; gqpt2 = sqrt(3/5); gqpt3 = -sqrt(3/5); gqpt = [gqpt1,gqpt2,gqpt3];
        gqwt1 = 8/9; gqwt2 = 5/9; gqwt3 = 5/9; gqwt = [gqwt1,gqwt2,gqwt3];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:size(elements,1) % j is the element index
    
    EleCrackTipOrNot = 0; EleCrackHOrNot = 0;
    StrainWithinEachElementGausspoint = zeros(4,4);
    
    % ----- Find four corner points -------
    point1x = coordinates(elements(j,1),1); point1y = coordinates(elements(j,1),2);
    point2x = coordinates(elements(j,2),1); point2y = coordinates(elements(j,2),2);
    point3x = coordinates(elements(j,3),1); point3y = coordinates(elements(j,3),2);
    point4x = coordinates(elements(j,4),1); point4y = coordinates(elements(j,4),2);
    
    % ------ Find mid points 5/6/7/8 -------
    if elements(j,5) ~= 0
        point5x = coordinates(elements(j,5),1); point5y = coordinates(elements(j,5),2);
    else, point5x = 0; point5y = 0; end
    if elements(j,6) ~= 0
        point6x = coordinates(elements(j,6),1); point6y = coordinates(elements(j,6),2);
    else, point6x = 0; point6y = 0; end
    if elements(j,7) ~= 0
        point7x = coordinates(elements(j,7),1); point7y = coordinates(elements(j,7),2);
    else, point7x = 0; point7y = 0; end
    if elements(j,8) ~= 0
        point8x = coordinates(elements(j,8),1); point8y = coordinates(elements(j,8),2);
    else, point8x = 0; point8y = 0; end

    % ------ Calculate ksi and eta --------
    lMatrix = [ point1x*point1y point1x point1y 1;
                point2x*point2y point2x point2y 1;
                point3x*point3y point3x point3y 1;
                point4x*point4y point4x point4y 1 ];

    % ------ Find the element nodal indices ------
    lb = [-1;1;1;-1]; l = linsolve(lMatrix,lb);
    mb = [-1;-1;1;1]; m = linsolve(lMatrix,mb);

    % ------ Find the element nodal indices ------
    temp = [2*elements(j,1)-1 2*elements(j,1) 2*elements(j,2)-1 2*elements(j,2) ...
            2*elements(j,3)-1 2*elements(j,3) 2*elements(j,4)-1 2*elements(j,4) ...
            2*elements(j,5)-1 2*elements(j,5) 2*elements(j,6)-1 2*elements(j,6) ...
            2*elements(j,7)-1 2*elements(j,7) 2*elements(j,8)-1 2*elements(j,8)];
        
    % We don't want temp <= 0, instead, put them to the end
    for tempi = 5:8
        if elements(j,tempi) == 0
            temp(2*tempi-1) = 2*(FEMSize+1)-1;    
            temp(2*tempi)   = 2*(FEMSize+1);
        end
    end
    
    % ------ Find the enriched functions nodal indices ------
    if CrackOrNot == 1 % Crack with and without crack tip
        NHIndexNode = (FEMSize+1)*ones(8,1); NTipIndexNode = (FEMSize+1)*ones(8,1);
        for tempi = 1:8
            if elements(j,tempi) ~= 0
                [tempindex,~] = find(EnrHAndTipEleIndex(:,1)==elements(j,tempi));
                if isempty(tempindex) ~= 1
                    NHIndexNode(tempi) = EnrHAndTipEleIndex(tempindex,2);
                    EleCrackHOrNot = EleCrackHOrNot + 1;
                end
                if CrackTipOrNot == 0
                [tempindexTip,~] = find(EnrTipEleIndex(:,1)==elements(j,tempi));
                    if isempty(tempindexTip) ~= 1
                        NTipIndexNode(tempi) = EnrHAndTipEleIndex(tempindexTip,2);
                        EleCrackTipOrNot = EleCrackTipOrNot + 1;
                    end
                end
            end
        end
        
        deltaPt5 = 0; deltaPt6 = 0; deltaPt7 = 0; deltaPt8 = 0;
        if elements(j,5) ~= 0; deltaPt5 = 1; end
        if elements(j,6) ~= 0; deltaPt6 = 1; end
        if elements(j,7) ~= 0; deltaPt7 = 1; end
        if elements(j,8) ~= 0; deltaPt8 = 1; end
        
        if EleCrackHOrNot > 0
            tempNH = [  2*NHIndexNode(1)-1 2*NHIndexNode(1) 2*NHIndexNode(2)-1 2*NHIndexNode(2) ...
                        2*NHIndexNode(3)-1 2*NHIndexNode(3) 2*NHIndexNode(4)-1 2*NHIndexNode(4) ...
                        2*NHIndexNode(5)-1 2*NHIndexNode(5) 2*NHIndexNode(6)-1 2*NHIndexNode(6) ...
                        2*NHIndexNode(7)-1 2*NHIndexNode(7) 2*NHIndexNode(8)-1 2*NHIndexNode(8)];
            temp = [temp tempNH];

            if EleCrackTipOrNot > 0  
                tempNTip = [2*NTipIndexNode(1)-1 2*NTipIndexNode(1) 2*(NTipIndexNode(1)+1)-1 2*(NTipIndexNode(1)+1) ...
                            2*(NTipIndexNode(1)+2)-1 2*(NTipIndexNode(1)+2) 2*(NTipIndexNode(1)+3)-1 2*(NTipIndexNode(1)+3) ...
                            2*NTipIndexNode(2)-1 2*NTipIndexNode(2) 2*(NTipIndexNode(2)+1)-1 2*(NTipIndexNode(2)+1) ...
                            2*(NTipIndexNode(2)+2)-1 2*(NTipIndexNode(2)+2) 2*(NTipIndexNode(2)+3)-1 2*(NTipIndexNode(2)+3) ...
                            2*NTipIndexNode(3)-1 2*NTipIndexNode(3) 2*(NTipIndexNode(3)+1)-1 2*(NTipIndexNode(3)+1) ...
                            2*(NTipIndexNode(3)+2)-1 2*(NTipIndexNode(3)+2) 2*(NTipIndexNode(3)+3)-1 2*(NTipIndexNode(3)+3) ...
                            2*NTipIndexNode(4)-1 2*NTipIndexNode(4) 2*(NTipIndexNode(4)+1)-1 2*(NTipIndexNode(4)+1) ...
                            2*(NTipIndexNode(4)+2)-1 2*(NTipIndexNode(4)+2) 2*(NTipIndexNode(4)+3)-1 2*(NTipIndexNode(4)+3) ...
                            2*NTipIndexNode(5)-1 2*NTipIndexNode(5) 2*(NTipIndexNode(5)+1)-1 2*(NTipIndexNode(5)+1) ...
                            2*(NTipIndexNode(5)+2)-1 2*(NTipIndexNode(5)+2) 2*(NTipIndexNode(5)+3)-1 2*(NTipIndexNode(5)+3) ...
                            2*NTipIndexNode(6)-1 2*NTipIndexNode(6)  2*(NTipIndexNode(6)+1)-1 2*(NTipIndexNode(6)+1) ...
                            2*(NTipIndexNode(6)+2)-1 2*(NTipIndexNode(6)+2)  2*(NTipIndexNode(6)+3)-1 2*(NTipIndexNode(6)+3) ...
                            2*NTipIndexNode(7)-1 2*NTipIndexNode(7) 2*(NTipIndexNode(7)+1)-1 2*(NTipIndexNode(7)+1) ...
                            2*(NTipIndexNode(7)+2)-1 2*(NTipIndexNode(7)+2)  2*(NTipIndexNode(7)+3)-1 2*(NTipIndexNode(7)+3) ...
                            2*NTipIndexNode(8)-1 2*NTipIndexNode(8) 2*(NTipIndexNode(8)+1)-1 2*(NTipIndexNode(8)+1) ...
                            2*(NTipIndexNode(8)+2)-1 2*(NTipIndexNode(8)+2) 2*(NTipIndexNode(8)+3)-1 2*(NTipIndexNode(8)+3) ];
                temp = [temp tempNTip];
            end
        end
         
                
    end
    
    % ------ Set Gauss points ------
    pt1 = elements(j,1); pt2 = elements(j,3);
    pointOfx = zeros(length(gqwt),1); pointOfy = zeros(length(gqwt),1);
    for tempi = 1:length(gqwt)
       pointOfx(tempi) = gqpt(tempi)*0.5*(coordinates(pt2,1)-coordinates(pt1,1))+0.5*(coordinates(pt2,1)+coordinates(pt1,1));
       pointOfy(tempi) = gqpt(tempi)*0.5*(coordinates(pt2,2)-coordinates(pt1,2))+0.5*(coordinates(pt2,2)+coordinates(pt1,2));
    end
    for tempi = 1:length(pointOfx)
       for tempj = 1:length(pointOfy)
    % ------ Start four Gauss points ------
    %for tempi = 1:2
    %    for tempj = 1:2
            
            % ------ Calculate ksi and eta ------
            ksi = l(1)*pointOfx(tempi)*pointOfy(tempj) + l(2)*pointOfx(tempi) + l(3)*pointOfy(tempj) + l(4) ;
            eta = m(1)*pointOfx(tempi)*pointOfy(tempj) + m(2)*pointOfx(tempi) + m(3)*pointOfy(tempj) + m(4) ;
            % ksi = 2*tempj-3; eta = 2*tempi-3;
            % if (tempi == 1) && (tempj == 1)
            %     ksi = -1/sqrt(3); eta = -1/sqrt(3);
            % elseif (tempi == 1) && (tempj == 2)
            %     ksi = 1/sqrt(3); eta = -1/sqrt(3);
            % elseif (tempi == 2) && (tempj == 1)
            %     ksi = 1/sqrt(3); eta = 1/sqrt(3);
            % elseif (tempi == 2) && (tempj == 2)
            %     ksi = -1/sqrt(3); eta = 1/sqrt(3);
            % end
            
            % ------ Calculate N ------
            deltaPt5 = 0; deltaPt6 = 0; deltaPt7 = 0; deltaPt8 = 0;
            if elements(j,5) ~= 0; deltaPt5 = 1; end
            if elements(j,6) ~= 0; deltaPt6 = 1; end
            if elements(j,7) ~= 0; deltaPt7 = 1; end
            if elements(j,8) ~= 0; deltaPt8 = 1; end

            N5 = deltaPt5*0.5*(1+ksi)*(1-abs(eta));
            N6 = deltaPt6*0.5*(1+eta)*(1-abs(ksi));
            N7 = deltaPt7*0.5*(1-ksi)*(1-abs(eta));
            N8 = deltaPt8*0.5*(1-eta)*(1-abs(ksi));

            N1 = (1-ksi)*(1-eta)*0.25 - 0.5*(N7+N8);
            N2 = (1+ksi)*(1-eta)*0.25 - 0.5*(N8+N5);
            N3 = (1+ksi)*(1+eta)*0.25 - 0.5*(N5+N6);
            N4 = (1-ksi)*(1+eta)*0.25 - 0.5*(N6+N7);

            %%%%%%%%% Consider crack enriched basis %%%%%%%%%%%%%%%
            if CrackOrNot == 0 || EleCrackHOrNot == 0 % No crack
                H1 = 0; H2 = 0; H3 = 0; H4 = 0; H5 = 0; H6 = 0; H7 = 0; H8 = 0; 
                C1 = 0; C2 = 0; C3 = 0; C4 = 0;
            elseif CrackOrNot == 1 && EleCrackHOrNot > 0 % Crack both with and without crack tip
                % Compute NH function: NH = N*H; so we need to compute
                % H function, to compute H function, we need to compare
                % coordinates of Gauss points with other nodes
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
                H1 = 0; if temp1 * tempNH < 0; H1 = sign(temp1-tempNH); end 
                H2 = 0; if temp2 * tempNH < 0; H2 = sign(temp2-tempNH); end 
                H3 = 0; if temp3 * tempNH < 0; H3 = sign(temp3-tempNH); end
                H4 = 0; if temp4 * tempNH < 0; H4 = sign(temp4-tempNH); end                     
                H5 = 0; if temp5 * tempNH < 0; H5 = sign(temp5-tempNH); end
                H6 = 0; if temp6 * tempNH < 0; H6 = sign(temp6-tempNH); end                    
                H7 = 0; if temp7 * tempNH < 0; H7 = sign(temp7-tempNH); end
                H8 = 0; if temp8 * tempNH < 0; H8 = sign(temp8-tempNH); end  
                % [tempNH, temp1, temp2, temp3, temp4, CrackTip(1), pointOfx(tempi), H1, H2, H3, H4 ]
                 
                C1 = 0; C2 = 0; C3 = 0; C4 = 0;
                if EleCrackTipOrNot > 0 % Crack with crack tip
                    % Compute C11-C44 function: Cij = Ni*Cj;
                    tempr = sqrt( (pointOfx(tempi)-CrackTip(1))^2 + (pointOfy(tempj)-CrackTip(2))^2 );
                    temptheta = atan2( (pointOfx(tempi)-CrackTip(1)) , (pointOfy(tempj)-CrackTip(2)) );  % atan2(Y,X)
                    C1 = sqrt(tempr)*sin(0.5*temptheta); 
                    C2 = sqrt(tempr)*cos(0.5*temptheta);
                    C3 = sqrt(tempr)*sin(0.5*temptheta)*sin(temptheta); 
                    C4 = sqrt(tempr)*cos(0.5*temptheta)*cos(temptheta);
                    
                    deltaCNode1 = 0; deltaCNode2 = 0; deltaCNode3 = 0; deltaCNode4 = 0; deltaCNode5 = 0; deltaCNode6 = 0; deltaCNode7 = 0; deltaCNode8 = 0;
                    if NTipIndexNode(1) ~= (FEMSize+1); deltaCNode1 = 1; end
                    if NTipIndexNode(2) ~= (FEMSize+1); deltaCNode2 = 1; end
                    if NTipIndexNode(3) ~= (FEMSize+1); deltaCNode3 = 1; end
                    if NTipIndexNode(4) ~= (FEMSize+1); deltaCNode4 = 1; end
                    if NTipIndexNode(5) ~= (FEMSize+1); deltaCNode5 = 1; end
                    if NTipIndexNode(6) ~= (FEMSize+1); deltaCNode6 = 1; end
                    if NTipIndexNode(7) ~= (FEMSize+1); deltaCNode7 = 1; end
                    if NTipIndexNode(8) ~= (FEMSize+1); deltaCNode8 = 1; end
                
                end
            end
            
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
            if CrackOrNot == 0 || EleCrackHOrNot == 0 % No crack, don't change previous codes
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
            elseif CrackOrNot == 1 && EleCrackHOrNot > 0 && EleCrackTipOrNot > 0 % Crack but without crack tip
                   DN = [InvJ zeros(2,2); zeros(2,2) InvJ] * ...
                    [funDN1Dksi(ksi,eta,deltaPt7,deltaPt8) 0 funDN2Dksi(ksi,eta,deltaPt8,deltaPt5) 0 ...
                    funDN3Dksi(ksi,eta,deltaPt5,deltaPt6) 0 funDN4Dksi(ksi,eta,deltaPt6,deltaPt7) 0 ...
                    funDN5Dksi(ksi,eta,deltaPt5) 0 funDN6Dksi(ksi,eta,deltaPt6) 0 ...
                    funDN7Dksi(ksi,eta,deltaPt7) 0 funDN8Dksi(ksi,eta,deltaPt8) 0 ...
                    H1*funDN1Dksi(ksi,eta,deltaPt7,deltaPt8) 0 H2*funDN2Dksi(ksi,eta,deltaPt8,deltaPt5) 0 ...
                    H3*funDN3Dksi(ksi,eta,deltaPt5,deltaPt6) 0 H4*funDN4Dksi(ksi,eta,deltaPt6,deltaPt7) 0 ...
                    H5*funDN5Dksi(ksi,eta,deltaPt5) 0 H6*funDN6Dksi(ksi,eta,deltaPt6) 0 ...
                    H7*funDN7Dksi(ksi,eta,deltaPt7) 0 H8*funDN8Dksi(ksi,eta,deltaPt8) 0 ...
                    C1*funDN1Dksi(ksi,eta,deltaPt7,deltaPt8)+N1*funDC1Dksi(ksi,eta,CrackTip) 0  C2*funDN1Dksi(ksi,eta,deltaPt7,deltaPt8)+N1*funDC2Dksi(ksi,eta,CrackTip) 0 ...
                    C3*funDN1Dksi(ksi,eta,deltaPt7,deltaPt8)+N1*funDC3Dksi(ksi,eta,CrackTip) 0  C4*funDN1Dksi(ksi,eta,deltaPt7,deltaPt8)+N1*funDC4Dksi(ksi,eta,CrackTip) 0 ...
                    C1*funDN2Dksi(ksi,eta,deltaPt8,deltaPt5)+N2*funDC1Dksi(ksi,eta,CrackTip) 0  C2*funDN2Dksi(ksi,eta,deltaPt8,deltaPt5)+N2*funDC2Dksi(ksi,eta,CrackTip) 0 ...
                    C3*funDN2Dksi(ksi,eta,deltaPt8,deltaPt5)+N2*funDC3Dksi(ksi,eta,CrackTip) 0  C4*funDN2Dksi(ksi,eta,deltaPt8,deltaPt5)+N2*funDC4Dksi(ksi,eta,CrackTip) 0 ...
                    C1*funDN3Dksi(ksi,eta,deltaPt5,deltaPt6)+N3*funDC1Dksi(ksi,eta,CrackTip) 0  C2*funDN3Dksi(ksi,eta,deltaPt5,deltaPt6)+N3*funDC2Dksi(ksi,eta,CrackTip) 0 ...
                    C3*funDN3Dksi(ksi,eta,deltaPt5,deltaPt6)+N3*funDC3Dksi(ksi,eta,CrackTip) 0  C4*funDN3Dksi(ksi,eta,deltaPt5,deltaPt6)+N3*funDC4Dksi(ksi,eta,CrackTip) 0 ...
                    C1*funDN4Dksi(ksi,eta,deltaPt6,deltaPt7)+N4*funDC1Dksi(ksi,eta,CrackTip) 0  C2*funDN4Dksi(ksi,eta,deltaPt6,deltaPt7)+N4*funDC2Dksi(ksi,eta,CrackTip) 0 ...
                    C3*funDN4Dksi(ksi,eta,deltaPt6,deltaPt7)+N4*funDC3Dksi(ksi,eta,CrackTip) 0  C4*funDN4Dksi(ksi,eta,deltaPt6,deltaPt7)+N4*funDC4Dksi(ksi,eta,CrackTip) 0 ...
                    C1*funDN5Dksi(ksi,eta,deltaPt5)+N5*funDC1Dksi(ksi,eta,CrackTip) 0 C2*funDN5Dksi(ksi,eta,deltaPt5)+N5*funDC2Dksi(ksi,eta,CrackTip) 0 ...
                    C3*funDN5Dksi(ksi,eta,deltaPt5)+N5*funDC3Dksi(ksi,eta,CrackTip) 0 C4*funDN5Dksi(ksi,eta,deltaPt5)+N5*funDC4Dksi(ksi,eta,CrackTip) 0 ...
                    C1*funDN6Dksi(ksi,eta,deltaPt6)+N6*funDC1Dksi(ksi,eta,CrackTip) 0 C2*funDN6Dksi(ksi,eta,deltaPt6)+N6*funDC2Dksi(ksi,eta,CrackTip) 0 ...
                    C3*funDN6Dksi(ksi,eta,deltaPt6)+N6*funDC3Dksi(ksi,eta,CrackTip) 0 C4*funDN6Dksi(ksi,eta,deltaPt6)+N6*funDC4Dksi(ksi,eta,CrackTip) 0 ...
                    C1*funDN7Dksi(ksi,eta,deltaPt7)+N7*funDC1Dksi(ksi,eta,CrackTip) 0 C2*funDN7Dksi(ksi,eta,deltaPt7)+N7*funDC2Dksi(ksi,eta,CrackTip) 0 ...
                    C3*funDN7Dksi(ksi,eta,deltaPt7)+N7*funDC3Dksi(ksi,eta,CrackTip) 0 C4*funDN7Dksi(ksi,eta,deltaPt7)+N7*funDC4Dksi(ksi,eta,CrackTip) 0 ...
                    C1*funDN8Dksi(ksi,eta,deltaPt8)+N8*funDC1Dksi(ksi,eta,CrackTip) 0 C2*funDN8Dksi(ksi,eta,deltaPt8)+N8*funDC2Dksi(ksi,eta,CrackTip) 0 ...
                    C3*funDN8Dksi(ksi,eta,deltaPt8)+N8*funDC3Dksi(ksi,eta,CrackTip) 0 C4*funDN8Dksi(ksi,eta,deltaPt8)+N8*funDC4Dksi(ksi,eta,CrackTip) 0 ;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    funDN1Deta(ksi,eta,deltaPt7,deltaPt8) 0 funDN2Deta(ksi,eta,deltaPt8,deltaPt5) 0 ...
                    funDN3Deta(ksi,eta,deltaPt5,deltaPt6) 0 funDN4Deta(ksi,eta,deltaPt6,deltaPt7) 0 ...
                    funDN5Deta(ksi,eta,deltaPt5) 0 funDN6Deta(ksi,eta,deltaPt6) 0 ...
                    funDN7Deta(ksi,eta,deltaPt7) 0 funDN8Deta(ksi,eta,deltaPt8) 0 ...
                    H1*funDN1Deta(ksi,eta,deltaPt7,deltaPt8) 0 H2*funDN2Deta(ksi,eta,deltaPt8,deltaPt5)   0 ...
                    H3*funDN3Deta(ksi,eta,deltaPt5,deltaPt6) 0 H4*funDN4Deta(ksi,eta,deltaPt6,deltaPt7)   0 ...
                    H5*funDN5Deta(ksi,eta,deltaPt5) 0 H6*funDN6Deta(ksi,eta,deltaPt6) 0 ...
                    H7*funDN7Deta(ksi,eta,deltaPt7) 0 H8*funDN8Deta(ksi,eta,deltaPt8) 0 ...
                    C1*funDN1Deta(ksi,eta,deltaPt7,deltaPt8)+N1*funDC1Deta(ksi,eta,CrackTip) 0  C2*funDN1Deta(ksi,eta,deltaPt7,deltaPt8)+N1*funDC2Deta(ksi,eta,CrackTip) 0 ...
                    C3*funDN1Deta(ksi,eta,deltaPt7,deltaPt8)+N1*funDC3Deta(ksi,eta,CrackTip) 0  C4*funDN1Deta(ksi,eta,deltaPt7,deltaPt8)+N1*funDC4Deta(ksi,eta,CrackTip) 0 ...
                    C1*funDN2Deta(ksi,eta,deltaPt8,deltaPt5)+N2*funDC1Deta(ksi,eta,CrackTip) 0  C2*funDN2Deta(ksi,eta,deltaPt8,deltaPt5)+N2*funDC2Deta(ksi,eta,CrackTip) 0 ...
                    C3*funDN2Deta(ksi,eta,deltaPt8,deltaPt5)+N2*funDC3Deta(ksi,eta,CrackTip) 0  C4*funDN2Deta(ksi,eta,deltaPt8,deltaPt5)+N2*funDC4Deta(ksi,eta,CrackTip) 0 ...
                    C1*funDN3Deta(ksi,eta,deltaPt5,deltaPt6)+N3*funDC1Deta(ksi,eta,CrackTip) 0  C2*funDN3Deta(ksi,eta,deltaPt5,deltaPt6)+N3*funDC2Deta(ksi,eta,CrackTip) 0 ...
                    C3*funDN3Deta(ksi,eta,deltaPt5,deltaPt6)+N3*funDC3Deta(ksi,eta,CrackTip) 0  C4*funDN3Deta(ksi,eta,deltaPt5,deltaPt6)+N3*funDC4Deta(ksi,eta,CrackTip) 0 ...
                    C1*funDN4Deta(ksi,eta,deltaPt6,deltaPt7)+N4*funDC1Deta(ksi,eta,CrackTip) 0  C2*funDN4Deta(ksi,eta,deltaPt6,deltaPt7)+N4*funDC2Deta(ksi,eta,CrackTip) 0 ...
                    C3*funDN4Deta(ksi,eta,deltaPt6,deltaPt7)+N4*funDC3Deta(ksi,eta,CrackTip) 0  C4*funDN4Deta(ksi,eta,deltaPt6,deltaPt7)+N4*funDC4Deta(ksi,eta,CrackTip) 0 ...
                    C1*funDN5Deta(ksi,eta,deltaPt5)+N5*funDC1Deta(ksi,eta,CrackTip) 0 C2*funDN5Deta(ksi,eta,deltaPt5)+N5*funDC2Deta(ksi,eta,CrackTip) 0 ...
                    C3*funDN5Deta(ksi,eta,deltaPt5)+N5*funDC3Deta(ksi,eta,CrackTip) 0 C4*funDN5Deta(ksi,eta,deltaPt5)+N5*funDC4Deta(ksi,eta,CrackTip) 0 ...
                    C1*funDN6Deta(ksi,eta,deltaPt6)+N6*funDC1Deta(ksi,eta,CrackTip) 0 C2*funDN6Deta(ksi,eta,deltaPt6)+N6*funDC2Deta(ksi,eta,CrackTip) 0 ...
                    C3*funDN6Deta(ksi,eta,deltaPt6)+N6*funDC3Deta(ksi,eta,CrackTip) 0 C4*funDN6Deta(ksi,eta,deltaPt6)+N6*funDC4Deta(ksi,eta,CrackTip) 0 ...
                    C1*funDN7Deta(ksi,eta,deltaPt7)+N7*funDC1Deta(ksi,eta,CrackTip) 0 C2*funDN7Deta(ksi,eta,deltaPt7)+N7*funDC2Deta(ksi,eta,CrackTip) 0 ...
                    C3*funDN7Deta(ksi,eta,deltaPt7)+N7*funDC3Deta(ksi,eta,CrackTip) 0 C4*funDN7Deta(ksi,eta,deltaPt7)+N7*funDC4Deta(ksi,eta,CrackTip) 0 ...
                    C1*funDN8Deta(ksi,eta,deltaPt8)+N8*funDC1Deta(ksi,eta,CrackTip) 0 C2*funDN8Deta(ksi,eta,deltaPt8)+N8*funDC2Deta(ksi,eta,CrackTip) 0 ...
                    C3*funDN8Deta(ksi,eta,deltaPt8)+N8*funDC3Deta(ksi,eta,CrackTip) 0 C4*funDN8Deta(ksi,eta,deltaPt8)+N8*funDC4Deta(ksi,eta,CrackTip) 0 ;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    0 funDN1Dksi(ksi,eta,deltaPt7,deltaPt8) 0 funDN2Dksi(ksi,eta,deltaPt8,deltaPt5) ...
                    0 funDN3Dksi(ksi,eta,deltaPt5,deltaPt6) 0 funDN4Dksi(ksi,eta,deltaPt6,deltaPt7) ...
                    0 funDN5Dksi(ksi,eta,deltaPt5) 0 funDN6Dksi(ksi,eta,deltaPt6) ...
                    0 funDN7Dksi(ksi,eta,deltaPt7) 0 funDN8Dksi(ksi,eta,deltaPt8) ...
                    0 H1*funDN1Dksi(ksi,eta,deltaPt7,deltaPt8) 0 H2*funDN2Dksi(ksi,eta,deltaPt8,deltaPt5) ...
                    0 H3*funDN3Dksi(ksi,eta,deltaPt5,deltaPt6) 0 H4*funDN4Dksi(ksi,eta,deltaPt6,deltaPt7) ...
                    0 H5*funDN5Dksi(ksi,eta,deltaPt5) 0 H6*funDN6Dksi(ksi,eta,deltaPt6) ...
                    0 H7*funDN7Dksi(ksi,eta,deltaPt7) 0 H8*funDN8Dksi(ksi,eta,deltaPt8) ...
                    0 C1*funDN1Dksi(ksi,eta,deltaPt7,deltaPt8)+N1*funDC1Dksi(ksi,eta,CrackTip) 0  C2*funDN1Dksi(ksi,eta,deltaPt7,deltaPt8)+N1*funDC2Dksi(ksi,eta,CrackTip)  ...
                    0 C3*funDN1Dksi(ksi,eta,deltaPt7,deltaPt8)+N1*funDC3Dksi(ksi,eta,CrackTip) 0  C4*funDN1Dksi(ksi,eta,deltaPt7,deltaPt8)+N1*funDC4Dksi(ksi,eta,CrackTip)  ...
                    0 C1*funDN2Dksi(ksi,eta,deltaPt8,deltaPt5)+N2*funDC1Dksi(ksi,eta,CrackTip) 0  C2*funDN2Dksi(ksi,eta,deltaPt8,deltaPt5)+N2*funDC2Dksi(ksi,eta,CrackTip)  ...
                    0 C3*funDN2Dksi(ksi,eta,deltaPt8,deltaPt5)+N2*funDC3Dksi(ksi,eta,CrackTip) 0  C4*funDN2Dksi(ksi,eta,deltaPt8,deltaPt5)+N2*funDC4Dksi(ksi,eta,CrackTip)  ...
                    0 C1*funDN3Dksi(ksi,eta,deltaPt5,deltaPt6)+N3*funDC1Dksi(ksi,eta,CrackTip) 0  C2*funDN3Dksi(ksi,eta,deltaPt5,deltaPt6)+N3*funDC2Dksi(ksi,eta,CrackTip)  ...
                    0 C3*funDN3Dksi(ksi,eta,deltaPt5,deltaPt6)+N3*funDC3Dksi(ksi,eta,CrackTip) 0  C4*funDN3Dksi(ksi,eta,deltaPt5,deltaPt6)+N3*funDC4Dksi(ksi,eta,CrackTip)  ...
                    0 C1*funDN4Dksi(ksi,eta,deltaPt6,deltaPt7)+N4*funDC1Dksi(ksi,eta,CrackTip) 0  C2*funDN4Dksi(ksi,eta,deltaPt6,deltaPt7)+N4*funDC2Dksi(ksi,eta,CrackTip)  ...
                    0 C3*funDN4Dksi(ksi,eta,deltaPt6,deltaPt7)+N4*funDC3Dksi(ksi,eta,CrackTip) 0  C4*funDN4Dksi(ksi,eta,deltaPt6,deltaPt7)+N4*funDC4Dksi(ksi,eta,CrackTip)  ...
                    0 C1*funDN5Dksi(ksi,eta,deltaPt5)+N5*funDC1Dksi(ksi,eta,CrackTip) 0 C2*funDN5Dksi(ksi,eta,deltaPt5)+N5*funDC2Dksi(ksi,eta,CrackTip)   ...
                    0 C3*funDN5Dksi(ksi,eta,deltaPt5)+N5*funDC3Dksi(ksi,eta,CrackTip) 0 C4*funDN5Dksi(ksi,eta,deltaPt5)+N5*funDC4Dksi(ksi,eta,CrackTip)   ...
                    0 C1*funDN6Dksi(ksi,eta,deltaPt6)+N6*funDC1Dksi(ksi,eta,CrackTip) 0 C2*funDN6Dksi(ksi,eta,deltaPt6)+N6*funDC2Dksi(ksi,eta,CrackTip)   ...
                    0 C3*funDN6Dksi(ksi,eta,deltaPt6)+N6*funDC3Dksi(ksi,eta,CrackTip) 0 C4*funDN6Dksi(ksi,eta,deltaPt6)+N6*funDC4Dksi(ksi,eta,CrackTip)   ...
                    0 C1*funDN7Dksi(ksi,eta,deltaPt7)+N7*funDC1Dksi(ksi,eta,CrackTip) 0 C2*funDN7Dksi(ksi,eta,deltaPt7)+N7*funDC2Dksi(ksi,eta,CrackTip)   ...
                    0 C3*funDN7Dksi(ksi,eta,deltaPt7)+N7*funDC3Dksi(ksi,eta,CrackTip) 0 C4*funDN7Dksi(ksi,eta,deltaPt7)+N7*funDC4Dksi(ksi,eta,CrackTip)   ...
                    0 C1*funDN8Dksi(ksi,eta,deltaPt8)+N8*funDC1Dksi(ksi,eta,CrackTip) 0 C2*funDN8Dksi(ksi,eta,deltaPt8)+N8*funDC2Dksi(ksi,eta,CrackTip)   ...
                    0 C3*funDN8Dksi(ksi,eta,deltaPt8)+N8*funDC3Dksi(ksi,eta,CrackTip) 0 C4*funDN8Dksi(ksi,eta,deltaPt8)+N8*funDC4Dksi(ksi,eta,CrackTip)  ;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    0 funDN1Deta(ksi,eta,deltaPt7,deltaPt8) 0 funDN2Deta(ksi,eta,deltaPt8,deltaPt5) ...
                    0 funDN3Deta(ksi,eta,deltaPt5,deltaPt6) 0 funDN4Deta(ksi,eta,deltaPt6,deltaPt7) ...
                    0 funDN5Deta(ksi,eta,deltaPt5) 0 funDN6Deta(ksi,eta,deltaPt6) ...
                    0 funDN7Deta(ksi,eta,deltaPt7) 0 funDN8Deta(ksi,eta,deltaPt8) ...
                    0 H1*funDN1Deta(ksi,eta,deltaPt7,deltaPt8) 0 H2*funDN2Deta(ksi,eta,deltaPt8,deltaPt5) ...
                    0 H3*funDN3Deta(ksi,eta,deltaPt5,deltaPt6) 0 H4*funDN4Deta(ksi,eta,deltaPt6,deltaPt7) ...
                    0 H5*funDN5Deta(ksi,eta,deltaPt5) 0 H6*funDN6Deta(ksi,eta,deltaPt6) ...
                    0 H7*funDN7Deta(ksi,eta,deltaPt7) 0 H8*funDN8Deta(ksi,eta,deltaPt8) ...
                    0 C1*funDN1Deta(ksi,eta,deltaPt7,deltaPt8)+N1*funDC1Deta(ksi,eta,CrackTip) 0  C2*funDN1Deta(ksi,eta,deltaPt7,deltaPt8)+N1*funDC2Deta(ksi,eta,CrackTip)  ...
                    0 C3*funDN1Deta(ksi,eta,deltaPt7,deltaPt8)+N1*funDC3Deta(ksi,eta,CrackTip) 0  C4*funDN1Deta(ksi,eta,deltaPt7,deltaPt8)+N1*funDC4Deta(ksi,eta,CrackTip)  ...
                    0 C1*funDN2Deta(ksi,eta,deltaPt8,deltaPt5)+N2*funDC1Deta(ksi,eta,CrackTip) 0  C2*funDN2Deta(ksi,eta,deltaPt8,deltaPt5)+N2*funDC2Deta(ksi,eta,CrackTip)  ...
                    0 C3*funDN2Deta(ksi,eta,deltaPt8,deltaPt5)+N2*funDC3Deta(ksi,eta,CrackTip) 0  C4*funDN2Deta(ksi,eta,deltaPt8,deltaPt5)+N2*funDC4Deta(ksi,eta,CrackTip)  ...
                    0 C1*funDN3Deta(ksi,eta,deltaPt5,deltaPt6)+N3*funDC1Deta(ksi,eta,CrackTip) 0  C2*funDN3Deta(ksi,eta,deltaPt5,deltaPt6)+N3*funDC2Deta(ksi,eta,CrackTip)  ...
                    0 C3*funDN3Deta(ksi,eta,deltaPt5,deltaPt6)+N3*funDC3Deta(ksi,eta,CrackTip) 0  C4*funDN3Deta(ksi,eta,deltaPt5,deltaPt6)+N3*funDC4Deta(ksi,eta,CrackTip)  ...
                    0 C1*funDN4Deta(ksi,eta,deltaPt6,deltaPt7)+N4*funDC1Deta(ksi,eta,CrackTip) 0  C2*funDN4Deta(ksi,eta,deltaPt6,deltaPt7)+N4*funDC2Deta(ksi,eta,CrackTip)  ...
                    0 C3*funDN4Deta(ksi,eta,deltaPt6,deltaPt7)+N4*funDC3Deta(ksi,eta,CrackTip) 0  C4*funDN4Deta(ksi,eta,deltaPt6,deltaPt7)+N4*funDC4Deta(ksi,eta,CrackTip)  ...
                    0 C1*funDN5Deta(ksi,eta,deltaPt5)+N5*funDC1Deta(ksi,eta,CrackTip) 0 C2*funDN5Deta(ksi,eta,deltaPt5)+N5*funDC2Deta(ksi,eta,CrackTip)   ...
                    0 C3*funDN5Deta(ksi,eta,deltaPt5)+N5*funDC3Deta(ksi,eta,CrackTip) 0 C4*funDN5Deta(ksi,eta,deltaPt5)+N5*funDC4Deta(ksi,eta,CrackTip)   ...
                    0 C1*funDN6Deta(ksi,eta,deltaPt6)+N6*funDC1Deta(ksi,eta,CrackTip) 0 C2*funDN6Deta(ksi,eta,deltaPt6)+N6*funDC2Deta(ksi,eta,CrackTip)   ...
                    0 C3*funDN6Deta(ksi,eta,deltaPt6)+N6*funDC3Deta(ksi,eta,CrackTip) 0 C4*funDN6Deta(ksi,eta,deltaPt6)+N6*funDC4Deta(ksi,eta,CrackTip)   ...
                    0 C1*funDN7Deta(ksi,eta,deltaPt7)+N7*funDC1Deta(ksi,eta,CrackTip) 0 C2*funDN7Deta(ksi,eta,deltaPt7)+N7*funDC2Deta(ksi,eta,CrackTip)   ...
                    0 C3*funDN7Deta(ksi,eta,deltaPt7)+N7*funDC3Deta(ksi,eta,CrackTip) 0 C4*funDN7Deta(ksi,eta,deltaPt7)+N7*funDC4Deta(ksi,eta,CrackTip)   ...
                    0 C1*funDN8Deta(ksi,eta,deltaPt8)+N8*funDC1Deta(ksi,eta,CrackTip) 0 C2*funDN8Deta(ksi,eta,deltaPt8)+N8*funDC2Deta(ksi,eta,CrackTip)   ...
                    0 C3*funDN8Deta(ksi,eta,deltaPt8)+N8*funDC3Deta(ksi,eta,CrackTip) 0 C4*funDN8Deta(ksi,eta,deltaPt8)+N8*funDC4Deta(ksi,eta,CrackTip) ];

            end
           
            StrainWithinEachElementGausspoint(length(gqpt)*(tempi-1)+tempj,1:4) = DN*U(temp);
            
        end
    end
    
    %%%%%%%%%%%%%%%%%%% Replace with 9-Gauss point %%%%%%%%%%%%%%%%%%%%%%
    switch GaussPtOrder
        case 3
            StrainGaussPt((length(gqpt)^2*(j-1)+1):length(gqpt)^2*j,1:4) = StrainWithinEachElementGausspoint;
            CoordGaussPt((length(gqpt)^2*(j-1)+1):length(gqpt)^2*j,1:2) = ...
                    [pointOfx(1),pointOfy(1);pointOfx(1),pointOfy(2);pointOfx(1),pointOfy(3);
                     pointOfx(2),pointOfy(1);pointOfx(2),pointOfy(2);pointOfx(2),pointOfy(3);
                     pointOfx(3),pointOfy(1);pointOfx(3),pointOfy(2);pointOfx(3),pointOfy(3)];
        %%%%%%%%%%%%%%%%%%% Comment following 4-Gauss point %%%%%%%%%%%%%%%%%%%%%%
        case 2
            StrainGaussPt((length(gqpt)^2*(j-1)+1):length(gqpt)^2*j,1:4) = StrainWithinEachElementGausspoint;
            CoordGaussPt((length(gqpt)^2*(j-1)+1):length(gqpt)^2*j,1:2) = ...
                    [pointOfx(1),pointOfy(1);pointOfx(1),pointOfy(2); 
                     pointOfx(2),pointOfy(1);pointOfx(2),pointOfy(2)];
                 
%             MatrixExtrapolation = [1+0.5*sqrt(3)  -0.5           1-0.5*sqrt(3)   -0.5;
%                                     -0.5           1+0.5*sqrt(3)  -0.5            1-0.5*sqrt(3);
%                                     1-0.5*sqrt(3)  -0.5           1+0.5*sqrt(3)   -0.5;
%                                     -0.5           1-0.5*sqrt(3)  -0.5            1+0.5*sqrt(3)];
%     
%             % ------ Nodal points strain extrapolation using Gauss points -----
%             StrainExxWithinEachElementNodalpoint =  MatrixExtrapolation * ...
%                 [StrainWithinEachElementGausspoint(1,1);
%                 StrainWithinEachElementGausspoint(2,1);
%                 StrainWithinEachElementGausspoint(3,1);
%                 StrainWithinEachElementGausspoint(4,1)];
% 
%             StrainExyWithinEachElementNodalpoint = MatrixExtrapolation * ...
%                 [StrainWithinEachElementGausspoint(1,2);
%                 StrainWithinEachElementGausspoint(2,2);
%                 StrainWithinEachElementGausspoint(3,2);
%                 StrainWithinEachElementGausspoint(4,2)];
% 
%             StrainEyxWithinEachElementNodalpoint = MatrixExtrapolation * ...
%                 [StrainWithinEachElementGausspoint(1,3);
%                 StrainWithinEachElementGausspoint(2,3);
%                 StrainWithinEachElementGausspoint(3,3);
%                 StrainWithinEachElementGausspoint(4,3)];
% 
%             StrainEyyWithinEachElementNodalpoint = MatrixExtrapolation * ...
%                 [StrainWithinEachElementGausspoint(1,4);
%                 StrainWithinEachElementGausspoint(2,4);
%                 StrainWithinEachElementGausspoint(3,4);
%                 StrainWithinEachElementGausspoint(4,4)];
% 
%             StrainWithinEachElementGausspoint(1,1) = StrainExxWithinEachElementNodalpoint(1);
%             StrainWithinEachElementGausspoint(2,1) = StrainExxWithinEachElementNodalpoint(2);
%             StrainWithinEachElementGausspoint(3,1) = StrainExxWithinEachElementNodalpoint(3);
%             StrainWithinEachElementGausspoint(4,1) = StrainExxWithinEachElementNodalpoint(4);
% 
%             StrainWithinEachElementGausspoint(1,2) = StrainExyWithinEachElementNodalpoint(1);
%             StrainWithinEachElementGausspoint(2,2) = StrainExyWithinEachElementNodalpoint(2);
%             StrainWithinEachElementGausspoint(3,2) = StrainExyWithinEachElementNodalpoint(3);
%             StrainWithinEachElementGausspoint(4,2) = StrainExyWithinEachElementNodalpoint(4);
% 
%             StrainWithinEachElementGausspoint(1,3) = StrainEyxWithinEachElementNodalpoint(1);
%             StrainWithinEachElementGausspoint(2,3) = StrainEyxWithinEachElementNodalpoint(2);
%             StrainWithinEachElementGausspoint(3,3) = StrainEyxWithinEachElementNodalpoint(3);
%             StrainWithinEachElementGausspoint(4,3) = StrainEyxWithinEachElementNodalpoint(4);
% 
%             StrainWithinEachElementGausspoint(1,4) = StrainEyyWithinEachElementNodalpoint(1);
%             StrainWithinEachElementGausspoint(2,4) = StrainEyyWithinEachElementNodalpoint(2);
%             StrainWithinEachElementGausspoint(3,4) = StrainEyyWithinEachElementNodalpoint(3);
%             StrainWithinEachElementGausspoint(4,4) = StrainEyyWithinEachElementNodalpoint(4);
% 
% 
%             % ------ Find the element nodal indices for strain ------
%             tempStrainIndex = [4*elements(j,1)-3  4*elements(j,1)-2  4*elements(j,1)-1  4*elements(j,1)  ...
%                                4*elements(j,2)-3  4*elements(j,2)-2  4*elements(j,2)-1  4*elements(j,2) ...
%                                4*elements(j,3)-3  4*elements(j,3)-2  4*elements(j,3)-1  4*elements(j,3) ...
%                                4*elements(j,4)-3  4*elements(j,4)-2  4*elements(j,4)-1  4*elements(j,4)];
% 
%             FStrain(tempStrainIndex) = FStrain(tempStrainIndex) + ...
%                [StrainWithinEachElementGausspoint(1,1); StrainWithinEachElementGausspoint(1,2); 
%                 StrainWithinEachElementGausspoint(1,3); StrainWithinEachElementGausspoint(1,4); 
%                 StrainWithinEachElementGausspoint(2,1); StrainWithinEachElementGausspoint(2,2); 
%                 StrainWithinEachElementGausspoint(2,3); StrainWithinEachElementGausspoint(2,4); 
%                 StrainWithinEachElementGausspoint(3,1); StrainWithinEachElementGausspoint(3,2); 
%                 StrainWithinEachElementGausspoint(3,3); StrainWithinEachElementGausspoint(3,4); 
%                 StrainWithinEachElementGausspoint(4,1); StrainWithinEachElementGausspoint(4,2); 
%                 StrainWithinEachElementGausspoint(4,3); StrainWithinEachElementGausspoint(4,4)];
% 
%             FStrainAvgTimes(tempStrainIndex) = FStrainAvgTimes(tempStrainIndex) + ones(16,1);
    end
    %%%%%%%%%%%%%%%%%%%%%%%% End of Comment %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

% switch GaussPtOrder
%     case 2
%         Strain = FStrain./FStrainAvgTimes;
%     case 3
%         Strain = zeros(size(coordinates,1)*4,1);
% end

% ====== Use gridfit and Gauss points to fit strain field ======
Coordxnodes = [min(coordinates(:,1)):winstepsize/(2^(LevelNo-1)):max(coordinates(:,1))]'; 
Coordynodes = [min(coordinates(:,2)):winstepsize/(2^(LevelNo-1)):max(coordinates(:,2))]';
Iblur_10 = gridfit(CoordGaussPt(:,1),CoordGaussPt(:,2),StrainGaussPt(:,1),Coordxnodes,Coordynodes); Iblur_10=Iblur_10';
Iblur_20 = gridfit(CoordGaussPt(:,1),CoordGaussPt(:,2),StrainGaussPt(:,2),Coordxnodes,Coordynodes); Iblur_20=Iblur_20';
Iblur_30 = gridfit(CoordGaussPt(:,1),CoordGaussPt(:,2),StrainGaussPt(:,3),Coordxnodes,Coordynodes); Iblur_30=Iblur_30';
Iblur_40 = gridfit(CoordGaussPt(:,1),CoordGaussPt(:,2),StrainGaussPt(:,4),Coordxnodes,Coordynodes); Iblur_40=Iblur_40';    
Strain = 0*[U(1:end-DIM*NodesNumPerEle);U(1:end-DIM*NodesNumPerEle)];             
for tempi = 1:size(coordinates,1)
        [row1,~] = find(Coordxnodes==coordinates(tempi,1));
        [row2,~] = find(Coordynodes==coordinates(tempi,2));
        Strain(4*tempi-3) = Iblur_10(row1,row2);
        Strain(4*tempi-2) = Iblur_20(row1,row2);
        Strain(4*tempi-1) = Iblur_30(row1,row2);
        Strain(4*tempi-0) = Iblur_40(row1,row2);
end
    

        
%% ========= subroutines for  FEM Q4 shape function derivatives ========

function DN5Dksi=funDN5Dksi(ksi,eta,deltaPt5)
DN5Dksi = deltaPt5*0.5*(1-abs(eta)) ;

function DN5Deta=funDN5Deta(ksi,eta,deltaPt5)
DN5Deta = deltaPt5*0.5*(1+ksi)*sign(-eta);

function DN6Dksi=funDN6Dksi(ksi,eta,deltaPt6)
DN6Dksi = deltaPt6*0.5*(1+eta)*sign(-ksi);

function DN6Deta=funDN6Deta(ksi,eta,deltaPt6)
DN6Deta = deltaPt6*0.5*(1-abs(ksi));

function DN7Dksi=funDN7Dksi(ksi,eta,deltaPt7)
DN7Dksi = deltaPt7*0.5*(-1)*(1-abs(eta));

function DN7Deta=funDN7Deta(ksi,eta,deltaPt7)
DN7Deta = deltaPt7*0.5*(1-ksi)*sign(-eta);

function DN8Dksi=funDN8Dksi(ksi,eta,deltaPt8)
DN8Dksi = deltaPt8*0.5*(1-eta)*sign(-ksi);

function DN8Deta=funDN8Deta(ksi,eta,deltaPt8)
DN8Deta = deltaPt8*0.5*(-1)*(1-abs(ksi));

function DN1Dksi = funDN1Dksi(ksi,eta,deltaPt7,deltaPt8)
DN1Dksi = -0.25*(1-eta)-0.5*((deltaPt7*0.5*(-1)*(1-abs(eta)))+deltaPt8*0.5*(1-eta)*sign(-ksi));

function DN1Deta = funDN1Deta(ksi,eta,deltaPt7,deltaPt8)
DN1Deta = -0.25*(1-ksi)-0.5*(deltaPt7*0.5*(1-ksi)*sign(-eta)+deltaPt8*0.5*(-1)*(1-abs(ksi)));

function DN2Dksi = funDN2Dksi(ksi,eta,deltaPt8,deltaPt5)
DN2Dksi = 0.25*(1-eta)-0.5*(deltaPt8*0.5*(1-eta)*sign(-ksi)+deltaPt5*0.5*(1-abs(eta)));

function DN2Deta = funDN2Deta(ksi,eta,deltaPt8,deltaPt5)
DN2Deta = -0.25*(1+ksi)-0.5*(deltaPt8*0.5*(-1)*(1-abs(ksi))+deltaPt5*0.5*(1+ksi)*sign(-eta));

function DN3Dksi = funDN3Dksi(ksi,eta,deltaPt5,deltaPt6)
DN3Dksi = 0.25*(1+eta)-0.5*(deltaPt5*0.5*(1-abs(eta))+deltaPt6*0.5*(1+eta)*sign(-ksi));

function DN3Deta = funDN3Deta(ksi,eta,deltaPt5,deltaPt6)
DN3Deta = 0.25*(1+ksi)-0.5*(deltaPt5*0.5*(1+ksi)*sign(-eta)+deltaPt6*0.5*(1-abs(ksi)));

function DN4Dksi = funDN4Dksi(ksi,eta,deltaPt6,deltaPt7)
DN4Dksi = -0.25*(1+eta)-0.5*(deltaPt6*0.5*(1+eta)*sign(-ksi)+deltaPt7*0.5*(-1)*(1-abs(eta)));

function DN4Deta = funDN4Deta(ksi,eta,deltaPt6,deltaPt7)
DN4Deta = 0.25*(1-ksi)-0.5*(deltaPt6*0.5*(1-abs(ksi))+deltaPt7*0.5*(1-ksi)*sign(-eta));

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



