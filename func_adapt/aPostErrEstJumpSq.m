function ErrEstJump = aPostErrEstJumpSq(coordinates,elements,j,eleGeneration,winsize,U)

ErrEstJump = 0; U = [U;0;0];

% Find neighboring elements for jump residuals
% We need function findEleNeighbors with +1, ==, -1 element generation;
% We also need to output the neighboring edges.

eleNeighborIndexAndEdge = findEleNeighbors(elements, coordinates, j ,eleGeneration, winsize);

% --------- Store format ---------
% for each row in the "eleNeighborIndexAndEdge"
% #neighboring element ID; #neighboring edge first node ID;  #neighboring edge second node ID;
% --------------------------------

% Gaussian quadrature parameter
gqpt1 = 0; gqpt2 = 0.538469; gqpt3 = -0.538469; gqpt4 = 0.90618; gqpt5 = -0.90618;
gqwt1 = 0.568889; gqwt2 = 0.478629; gqwt3 = 0.478629; gqwt4 = 0.236927; gqwt5 = 0.236927;

gqwt = [gqwt1,gqwt2,gqwt3,gqwt4,gqwt5];


%% ======== Left side =========

if eleNeighborIndexAndEdge(1,1) ~= 0
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ========= Element itself ================
    
    % ------ Find four corner points ------
    point1x = coordinates(elements(j,1),1);
    point1y = coordinates(elements(j,1),2);
    point2x = coordinates(elements(j,2),1);
    point2y = coordinates(elements(j,2),2);
    point3x = coordinates(elements(j,3),1);
    point3y = coordinates(elements(j,3),2);
    point4x = coordinates(elements(j,4),1);
    point4y = coordinates(elements(j,4),2);
    
    lengthOfElement = point3x-point1x;
    
    % ------ Find mid points 5/6/7/8 -------
    if elements(j,5) ~= 0
        point5x = coordinates(elements(j,5),1); point5y = coordinates(elements(j,5),2);
    else
        point5x = 0; point5y = 0;
    end
    
    if elements(j,6) ~= 0
        point6x = coordinates(elements(j,6),1); point6y = coordinates(elements(j,6),2);
    else
        point6x = 0; point6y = 0;
    end
    
    if elements(j,7) ~= 0
        point7x = coordinates(elements(j,7),1); point7y = coordinates(elements(j,7),2);
    else
        point7x = 0; point7y = 0;
    end
    
    if elements(j,8) ~= 0
        point8x = coordinates(elements(j,8),1); point8y = coordinates(elements(j,8),2);
    else
        point8x = 0; point8y = 0;
    end
    
    % ------ Calculate ksi and eta --------
    lMatrix = [point1x*point1y point1x point1y 1;
        point2x*point2y point2x point2y 1;
        point3x*point3y point3x point3y 1;
        point4x*point4y point4x point4y 1];
    
    lb = [-1;1;1;-1];
    l = linsolve(lMatrix,lb);
    
    mb = [-1;-1;1;1];
    m = linsolve(lMatrix,mb);
     
    % ------ Find the element nodal indices ------
    temp = [2*elements(j,1)-1 2*elements(j,1) 2*elements(j,2)-1 2*elements(j,2) ...
        2*elements(j,3)-1 2*elements(j,3) 2*elements(j,4)-1 2*elements(j,4) ...
        2*elements(j,5)-1 2*elements(j,5) 2*elements(j,6)-1 2*elements(j,6) ...
        2*elements(j,7)-1 2*elements(j,7) 2*elements(j,8)-1 2*elements(j,8)];
    
    % We don't want temp <= 0, instead, put them to the end
    for tempi = 5:8
        if elements(j,tempi) == 0
            temp(2*tempi-1) = 2*(size(coordinates,1)+1)-1;
            temp(2*tempi)   = 2*(size(coordinates,1)+1);
        end
    end
    
    pt1 = eleNeighborIndexAndEdge(1,2);
    pt2 = eleNeighborIndexAndEdge(1,3);
    pointOfx = zeros(5,1); pointOfy = zeros(5,1);
    
    pointOfx(1) = gqpt1*0.5*(coordinates(pt2,1)-coordinates(pt1,1))+0.5*(coordinates(pt2,1)+coordinates(pt1,1));
    pointOfy(1) = gqpt1*0.5*(coordinates(pt2,2)-coordinates(pt1,2))+0.5*(coordinates(pt2,2)+coordinates(pt1,2));
    % pointOfx(2) = gqpt2*0.5*(coordinates(pt2,1)-coordinates(pt1,1))+0.5*(coordinates(pt2,1)+coordinates(pt1,1));
    pointOfy(2) = gqpt2*0.5*(coordinates(pt2,2)-coordinates(pt1,2))+0.5*(coordinates(pt2,2)+coordinates(pt1,2));
    % pointOfx(3) = gqpt3*0.5*(coordinates(pt2,1)-coordinates(pt1,1))+0.5*(coordinates(pt2,1)+coordinates(pt1,1));
    pointOfy(3) = gqpt3*0.5*(coordinates(pt2,2)-coordinates(pt1,2))+0.5*(coordinates(pt2,2)+coordinates(pt1,2));
    % pointOfx(4) = gqpt4*0.5*(coordinates(pt2,1)-coordinates(pt1,1))+0.5*(coordinates(pt2,1)+coordinates(pt1,1));
    pointOfy(4) = gqpt4*0.5*(coordinates(pt2,2)-coordinates(pt1,2))+0.5*(coordinates(pt2,2)+coordinates(pt1,2));
    % pointOfx(5) = gqpt5*0.5*(coordinates(pt2,1)-coordinates(pt1,1))+0.5*(coordinates(pt2,1)+coordinates(pt1,1));
    pointOfy(5) = gqpt5*0.5*(coordinates(pt2,2)-coordinates(pt1,2))+0.5*(coordinates(pt2,2)+coordinates(pt1,2));
    
    DUDX1 = zeros(size(pointOfx,1),size(pointOfy,1));
    DVDX1 = zeros(size(pointOfx,1),size(pointOfy,1));
    
    for tempi = 1 % :size(pointOfx,1)
        for tempj = 1: size(pointOfy,1)
            
            % ------ Calculate ksi and eta ------
            ksi = l(1)*pointOfx(tempi)*pointOfy(tempj) + l(2)*pointOfx(tempi) + l(3)*pointOfy(tempj) + l(4) ;
            eta = m(1)*pointOfx(tempi)*pointOfy(tempj) + m(2)*pointOfx(tempi) + m(3)*pointOfy(tempj) + m(4) ;
            
            % ------ Calculate N ------
            deltaPt5 = 0; deltaPt6 = 0; deltaPt7 = 0; deltaPt8 = 0;
            if elements(j,5) ~= 0
                deltaPt5 = 1;
            end
            if elements(j,6) ~= 0
                deltaPt6 = 1;
            end
            if elements(j,7) ~= 0
                deltaPt7 = 1;
            end
            if elements(j,8) ~= 0
                deltaPt8 = 1;
            end
            
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
            
            % ------- Calculate DU/Dx ---------
            tempDUDX = DN*U(temp);
            DUDX1(tempi,tempj) = tempDUDX(1);
            DVDX1(tempi,tempj) = tempDUDX(3);
            
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ========= Left side 1 element ===========
    
    elementtemp = eleNeighborIndexAndEdge(1,1);
    
    % ------ Find four corner points ------
    point1x = coordinates(elements(eleNeighborIndexAndEdge(1,1),1),1);
    point1y = coordinates(elements(eleNeighborIndexAndEdge(1,1),1),2);
    point2x = coordinates(elements(eleNeighborIndexAndEdge(1,1),2),1);
    point2y = coordinates(elements(eleNeighborIndexAndEdge(1,1),2),2);
    point3x = coordinates(elements(eleNeighborIndexAndEdge(1,1),3),1);
    point3y = coordinates(elements(eleNeighborIndexAndEdge(1,1),3),2);
    point4x = coordinates(elements(eleNeighborIndexAndEdge(1,1),4),1);
    point4y = coordinates(elements(eleNeighborIndexAndEdge(1,1),4),2);
    
    % ------ Find mid points 5/6/7/8 -------
    if elements(eleNeighborIndexAndEdge(1,1),5) ~= 0
        point5x = coordinates(elements(eleNeighborIndexAndEdge(1,1),5),1); 
        point5y = coordinates(elements(eleNeighborIndexAndEdge(1,1),5),2);
    else
        point5x = 0; point5y = 0;
    end
    
    if elements(eleNeighborIndexAndEdge(1,1),6) ~= 0
        point6x = coordinates(elements(eleNeighborIndexAndEdge(1,1),6),1); 
        point6y = coordinates(elements(eleNeighborIndexAndEdge(1,1),6),2);
    else
        point6x = 0; point6y = 0;
    end
    
    if elements(eleNeighborIndexAndEdge(1,1),7) ~= 0
        point7x = coordinates(elements(eleNeighborIndexAndEdge(1,1),7),1); 
        point7y = coordinates(elements(eleNeighborIndexAndEdge(1,1),7),2);
    else
        point7x = 0; point7y = 0;
    end
    
    if elements(eleNeighborIndexAndEdge(1,1),8) ~= 0
        point8x = coordinates(elements(eleNeighborIndexAndEdge(1,1),8),1); 
        point8y = coordinates(elements(eleNeighborIndexAndEdge(1,1),8),2);
    else
        point8x = 0; point8y = 0;
    end
    
    % ------ Calculate ksi and eta --------
    lMatrix = [point1x*point1y point1x point1y 1;
        point2x*point2y point2x point2y 1;
        point3x*point3y point3x point3y 1;
        point4x*point4y point4x point4y 1];
    
    lb = [-1;1;1;-1];
    l = linsolve(lMatrix,lb);
    
    mb = [-1;-1;1;1];
    m = linsolve(lMatrix,mb);
    
    % ------ Find the element nodal indices ------
    temp = [2*elements(eleNeighborIndexAndEdge(1,1),1)-1 2*elements(eleNeighborIndexAndEdge(1,1),1) ...
        2*elements(eleNeighborIndexAndEdge(1,1),2)-1 2*elements(eleNeighborIndexAndEdge(1,1),2) ...
        2*elements(eleNeighborIndexAndEdge(1,1),3)-1 2*elements(eleNeighborIndexAndEdge(1,1),3) ...
        2*elements(eleNeighborIndexAndEdge(1,1),4)-1 2*elements(eleNeighborIndexAndEdge(1,1),4) ...
        2*elements(eleNeighborIndexAndEdge(1,1),5)-1 2*elements(eleNeighborIndexAndEdge(1,1),5) ...
        2*elements(eleNeighborIndexAndEdge(1,1),6)-1 2*elements(eleNeighborIndexAndEdge(1,1),6) ...
        2*elements(eleNeighborIndexAndEdge(1,1),7)-1 2*elements(eleNeighborIndexAndEdge(1,1),7) ...
        2*elements(eleNeighborIndexAndEdge(1,1),8)-1 2*elements(eleNeighborIndexAndEdge(1,1),8)];
    
    % We don't want temp <= 0, instead, put them to the end
    for tempi = 5:8
        if elements(eleNeighborIndexAndEdge(1,1),tempi) == 0
            temp(2*tempi-1) = 2*(size(coordinates,1)+1)-1;
            temp(2*tempi)   = 2*(size(coordinates,1)+1);
        end
    end
    
    
    DUDXLeft1 = zeros(size(pointOfx,1),size(pointOfy,1));
    DVDXLeft1 = zeros(size(pointOfx,1),size(pointOfy,1));
    
    for tempi = 1  %:size(pointOfx,1)
        for tempj = 1:size(pointOfy,1)
            
            % ------ Calculate ksi and eta -------
            ksi = l(1)*pointOfx(tempi)*pointOfy(tempj) + l(2)*pointOfx(tempi) + l(3)*pointOfy(tempj) + l(4);
            eta = m(1)*pointOfx(tempi)*pointOfy(tempj) + m(2)*pointOfx(tempi) + m(3)*pointOfy(tempj) + m(4);
            
            % ------ Calculate N ------
            deltaPt5 = 0; deltaPt6 = 0; deltaPt7 = 0; deltaPt8 = 0;
            if elements(eleNeighborIndexAndEdge(1,1),5) ~= 0
                deltaPt5 = 1;
            end
            if elements(eleNeighborIndexAndEdge(1,1),6) ~= 0
                deltaPt6 = 1;
            end
            if elements(eleNeighborIndexAndEdge(1,1),7) ~= 0
                deltaPt7 = 1;
            end
            if elements(eleNeighborIndexAndEdge(1,1),8) ~= 0
                deltaPt8 = 1;
            end
            
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
            
            % ------- Calculate DU/Dx ---------
            tempDUDX = DN*U(temp);
            DUDXLeft1(tempi,tempj) = tempDUDX(1);
            DVDXLeft1(tempi,tempj) = tempDUDX(3);
            
        end
    end
    
end


% ======== Left side 2 neighbor =========
if eleNeighborIndexAndEdge(2,1) ~= 0
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ========= Element itself ================
    
    % ------ Find four corner points ------
    point1x = coordinates(elements(j,1),1);
    point1y = coordinates(elements(j,1),2);
    point2x = coordinates(elements(j,2),1);
    point2y = coordinates(elements(j,2),2);
    point3x = coordinates(elements(j,3),1);
    point3y = coordinates(elements(j,3),2);
    point4x = coordinates(elements(j,4),1);
    point4y = coordinates(elements(j,4),2);
    
    % ------ Find mid points 5/6/7/8 -------
    if elements(j,5) ~= 0
        point5x = coordinates(elements(j,5),1); point5y = coordinates(elements(j,5),2);
    else
        point5x = 0; point5y = 0;
    end
    
    if elements(j,6) ~= 0
        point6x = coordinates(elements(j,6),1); point6y = coordinates(elements(j,6),2);
    else
        point6x = 0; point6y = 0;
    end
    
    if elements(j,7) ~= 0
        point7x = coordinates(elements(j,7),1); point7y = coordinates(elements(j,7),2);
    else
        point7x = 0; point7y = 0;
    end
    
    if elements(j,8) ~= 0
        point8x = coordinates(elements(j,8),1); point8y = coordinates(elements(j,8),2);
    else
        point8x = 0; point8y = 0;
    end
    
    % ------ Calculate ksi and eta --------
    lMatrix = [point1x*point1y point1x point1y 1;
        point2x*point2y point2x point2y 1;
        point3x*point3y point3x point3y 1;
        point4x*point4y point4x point4y 1];
    
    lb = [-1;1;1;-1];
    l = linsolve(lMatrix,lb);
    
    mb = [-1;-1;1;1];
    m = linsolve(lMatrix,mb);
     
    % ------ Find the element nodal indices ------
    temp = [2*elements(j,1)-1 2*elements(j,1) 2*elements(j,2)-1 2*elements(j,2) ...
        2*elements(j,3)-1 2*elements(j,3) 2*elements(j,4)-1 2*elements(j,4) ...
        2*elements(j,5)-1 2*elements(j,5) 2*elements(j,6)-1 2*elements(j,6) ...
        2*elements(j,7)-1 2*elements(j,7) 2*elements(j,8)-1 2*elements(j,8)];
    
    % We don't want temp <= 0, instead, put them to the end
    for tempi = 5:8
        if elements(j,tempi) == 0
            temp(2*tempi-1) = 2*(size(coordinates,1)+1)-1;
            temp(2*tempi)   = 2*(size(coordinates,1)+1);
        end
    end
    
    pt1 = eleNeighborIndexAndEdge(2,2);
    pt2 = eleNeighborIndexAndEdge(2,3);
    pointOfx = zeros(5,1); pointOfy = zeros(5,1);
    
    % Following pointOfx(1:5) and pointOfy(1:5) are used for Gaussian Quadrature
    
    pointOfx(1) = gqpt1*0.5*(coordinates(pt2,1)-coordinates(pt1,1))+0.5*(coordinates(pt2,1)+coordinates(pt1,1));
    pointOfy(1) = gqpt1*0.5*(coordinates(pt2,2)-coordinates(pt1,2))+0.5*(coordinates(pt2,2)+coordinates(pt1,2));
    % pointOfx(2) = gqpt2*0.5*(coordinates(pt2,1)-coordinates(pt1,1))+0.5*(coordinates(pt2,1)+coordinates(pt1,1));
    pointOfy(2) = gqpt2*0.5*(coordinates(pt2,2)-coordinates(pt1,2))+0.5*(coordinates(pt2,2)+coordinates(pt1,2));
    % pointOfx(3) = gqpt3*0.5*(coordinates(pt2,1)-coordinates(pt1,1))+0.5*(coordinates(pt2,1)+coordinates(pt1,1));
    pointOfy(3) = gqpt3*0.5*(coordinates(pt2,2)-coordinates(pt1,2))+0.5*(coordinates(pt2,2)+coordinates(pt1,2));
    % pointOfx(4) = gqpt4*0.5*(coordinates(pt2,1)-coordinates(pt1,1))+0.5*(coordinates(pt2,1)+coordinates(pt1,1));
    pointOfy(4) = gqpt4*0.5*(coordinates(pt2,2)-coordinates(pt1,2))+0.5*(coordinates(pt2,2)+coordinates(pt1,2));
    % pointOfx(5) = gqpt5*0.5*(coordinates(pt2,1)-coordinates(pt1,1))+0.5*(coordinates(pt2,1)+coordinates(pt1,1));
    pointOfy(5) = gqpt5*0.5*(coordinates(pt2,2)-coordinates(pt1,2))+0.5*(coordinates(pt2,2)+coordinates(pt1,2));
    
    DUDX2 = zeros(size(pointOfx,1),size(pointOfy,1));
    DVDX2 = zeros(size(pointOfx,1),size(pointOfy,1));
    
    for tempi = 1 % :size(pointOfx,1)
        for tempj = 1: size(pointOfy,1)
            
            % ------ Calculate ksi and eta ------
            ksi = l(1)*pointOfx(tempi)*pointOfy(tempj) + l(2)*pointOfx(tempi) + l(3)*pointOfy(tempj) + l(4) ;
            eta = m(1)*pointOfx(tempi)*pointOfy(tempj) + m(2)*pointOfx(tempi) + m(3)*pointOfy(tempj) + m(4) ;
            
            % ------ Calculate N ------
            deltaPt5 = 0; deltaPt6 = 0; deltaPt7 = 0; deltaPt8 = 0;
            if elements(j,5) ~= 0
                deltaPt5 = 1;
            end
            if elements(j,6) ~= 0
                deltaPt6 = 1;
            end
            if elements(j,7) ~= 0
                deltaPt7 = 1;
            end
            if elements(j,8) ~= 0
                deltaPt8 = 1;
            end
            
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
            
            % ------- Calculate DU/Dx ---------
            tempDUDX = DN*U(temp);
            DUDX2(tempi,tempj) = tempDUDX(1);
            DVDX2(tempi,tempj) = tempDUDX(3);
            
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ========= Left side element ===========
    
    % ------ Find four corner points ------
    point1x = coordinates(elements(eleNeighborIndexAndEdge(2,1),1),1);
    point1y = coordinates(elements(eleNeighborIndexAndEdge(2,1),1),2);
    point2x = coordinates(elements(eleNeighborIndexAndEdge(2,1),2),1);
    point2y = coordinates(elements(eleNeighborIndexAndEdge(2,1),2),2);
    point3x = coordinates(elements(eleNeighborIndexAndEdge(2,1),3),1);
    point3y = coordinates(elements(eleNeighborIndexAndEdge(2,1),3),2);
    point4x = coordinates(elements(eleNeighborIndexAndEdge(2,1),4),1);
    point4y = coordinates(elements(eleNeighborIndexAndEdge(2,1),4),2);
    
    % ------ Find mid points 5/6/7/8 -------
    if elements(eleNeighborIndexAndEdge(2,1),5) ~= 0
        point5x = coordinates(elements(eleNeighborIndexAndEdge(2,1),5),1); 
        point5y = coordinates(elements(eleNeighborIndexAndEdge(2,1),5),2);
    else
        point5x = 0; point5y = 0;
    end
    
    if elements(eleNeighborIndexAndEdge(2,1),6) ~= 0
        point6x = coordinates(elements(eleNeighborIndexAndEdge(2,1),6),1); 
        point6y = coordinates(elements(eleNeighborIndexAndEdge(2,1),6),2);
    else
        point6x = 0; point6y = 0;
    end
    
    if elements(eleNeighborIndexAndEdge(2,1),7) ~= 0
        point7x = coordinates(elements(eleNeighborIndexAndEdge(2,1),7),1); 
        point7y = coordinates(elements(eleNeighborIndexAndEdge(2,1),7),2);
    else
        point7x = 0; point7y = 0;
    end
    
    if elements(eleNeighborIndexAndEdge(2,1),8) ~= 0
        point8x = coordinates(elements(eleNeighborIndexAndEdge(2,1),8),1); 
        point8y = coordinates(elements(eleNeighborIndexAndEdge(2,1),8),2);
    else
        point8x = 0; point8y = 0;
    end
    
    % ------ Calculate ksi and eta --------
    lMatrix = [ point1x*point1y point1x point1y 1;
                point2x*point2y point2x point2y 1;
                point3x*point3y point3x point3y 1;
                point4x*point4y point4x point4y 1 ];
        
    lb = [-1;1;1;-1];
    l = linsolve(lMatrix,lb);
    
    mb = [-1;-1;1;1];
    m = linsolve(lMatrix,mb);
    
    % ------ Find the element nodal indices ------
    temp = [2*elements(eleNeighborIndexAndEdge(2,1),1)-1 2*elements(eleNeighborIndexAndEdge(2,1),1) ...
        2*elements(eleNeighborIndexAndEdge(2,1),2)-1 2*elements(eleNeighborIndexAndEdge(2,1),2) ...
        2*elements(eleNeighborIndexAndEdge(2,1),3)-1 2*elements(eleNeighborIndexAndEdge(2,1),3) ...
        2*elements(eleNeighborIndexAndEdge(2,1),4)-1 2*elements(eleNeighborIndexAndEdge(2,1),4) ...
        2*elements(eleNeighborIndexAndEdge(2,1),5)-1 2*elements(eleNeighborIndexAndEdge(2,1),5) ...
        2*elements(eleNeighborIndexAndEdge(2,1),6)-1 2*elements(eleNeighborIndexAndEdge(2,1),6) ...
        2*elements(eleNeighborIndexAndEdge(2,1),7)-1 2*elements(eleNeighborIndexAndEdge(2,1),7) ...
        2*elements(eleNeighborIndexAndEdge(2,1),8)-1 2*elements(eleNeighborIndexAndEdge(2,1),8)];
    
    % We don't want temp <= 0, instead, put them to the end
    for tempi = 5:8
        if elements(eleNeighborIndexAndEdge(2,1),tempi) == 0
            temp(2*tempi-1) = 2*(size(coordinates,1)+1)-1;
            temp(2*tempi)   = 2*(size(coordinates,1)+1);
        end
    end
    
    
    DUDXLeft2 = zeros(size(pointOfx,1),size(pointOfy,1));
    DVDXLeft2 = zeros(size(pointOfx,1),size(pointOfy,1));
    
    for tempi = 1  %:size(pointOfx,1)
        for tempj = 1:size(pointOfy,1)
            
            % ------ Calculate ksi and eta -------
            ksi = l(1)*pointOfx(tempi)*pointOfy(tempj) + l(2)*pointOfx(tempi) + l(3)*pointOfy(tempj) + l(4);
            eta = m(1)*pointOfx(tempi)*pointOfy(tempj) + m(2)*pointOfx(tempi) + m(3)*pointOfy(tempj) + m(4);
            
            % ------ Calculate N ------
            deltaPt5 = 0; deltaPt6 = 0; deltaPt7 = 0; deltaPt8 = 0;
            if elements(eleNeighborIndexAndEdge(2,1),5) ~= 0
                deltaPt5 = 1;
            end
            if elements(eleNeighborIndexAndEdge(2,1),6) ~= 0
                deltaPt6 = 1;
            end
            if elements(eleNeighborIndexAndEdge(2,1),7) ~= 0
                deltaPt7 = 1;
            end
            if elements(eleNeighborIndexAndEdge(2,1),8) ~= 0
                deltaPt8 = 1;
            end
            
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
            
            % ------- Calculate DU/Dx ---------
            tempDUDX = DN*U(temp);
            DUDXLeft2(tempi,tempj) = tempDUDX(1);
            DVDXLeft2(tempi,tempj) = tempDUDX(3);
            
        end
    end
    
end


% ============ Compute jump error ===============

if eleNeighborIndexAndEdge(1,1) ~= 0
    if eleNeighborIndexAndEdge(2,1) ~= 0
        for tempi = 1:length(gqwt)
            for tempj = 1:length(gqwt)
                ErrEstJump = ErrEstJump + gqwt(tempi)*gqwt(tempj)*lengthOfElement*0.5*((DUDX1(tempi,tempj)-DUDXLeft1(tempi,tempj))^2 ...
                             + (DVDX1(tempi,tempj)-DVDXLeft1(tempi,tempj))^2) + ...
                             gqwt(tempj)*gqwt(tempi)*lengthOfElement*0.5*((DUDX2(tempi,tempj)-DUDXLeft2(tempi,tempj))^2 ...
                             + (DVDX2(tempi,tempj)-DVDXLeft2(tempi,tempj))^2 );
            end
        end
    else
        for tempi = 1:length(gqwt)
            for tempj = 1:length(gqwt)
                ErrEstJump = ErrEstJump + gqwt(tempi)*gqwt(tempj)*lengthOfElement*((DUDX1(tempi,tempj)-DUDXLeft1(tempi,tempj))^2 ...
                             + (DVDX1(tempi,tempj)-DVDXLeft1(tempi,tempj))^2);
            end
        end
    end
end





%% % ======== Right side =========

if eleNeighborIndexAndEdge(3,1) ~= 0
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ========= Element itself ================
    
    % ------ Find four corner points ------
    point1x = coordinates(elements(j,1),1);
    point1y = coordinates(elements(j,1),2);
    point2x = coordinates(elements(j,2),1);
    point2y = coordinates(elements(j,2),2);
    point3x = coordinates(elements(j,3),1);
    point3y = coordinates(elements(j,3),2);
    point4x = coordinates(elements(j,4),1);
    point4y = coordinates(elements(j,4),2);
    
    lengthOfElement = point3x-point1x;
    
    % ------ Find mid points 5/6/7/8 -------
    if elements(j,5) ~= 0
        point5x = coordinates(elements(j,5),1); point5y = coordinates(elements(j,5),2);
    else
        point5x = 0; point5y = 0;
    end
    
    if elements(j,6) ~= 0
        point6x = coordinates(elements(j,6),1); point6y = coordinates(elements(j,6),2);
    else
        point6x = 0; point6y = 0;
    end
    
    if elements(j,7) ~= 0
        point7x = coordinates(elements(j,7),1); point7y = coordinates(elements(j,7),2);
    else
        point7x = 0; point7y = 0;
    end
    
    if elements(j,8) ~= 0
        point8x = coordinates(elements(j,8),1); point8y = coordinates(elements(j,8),2);
    else
        point8x = 0; point8y = 0;
    end
    
    % ------ Calculate ksi and eta --------
    lMatrix = [point1x*point1y point1x point1y 1;
        point2x*point2y point2x point2y 1;
        point3x*point3y point3x point3y 1;
        point4x*point4y point4x point4y 1];
    
    lb = [-1;1;1;-1];
    l = linsolve(lMatrix,lb);
    
    mb = [-1;-1;1;1];
    m = linsolve(lMatrix,mb);
     
    % ------ Find the element nodal indices ------
    temp = [2*elements(j,1)-1 2*elements(j,1) 2*elements(j,2)-1 2*elements(j,2) ...
        2*elements(j,3)-1 2*elements(j,3) 2*elements(j,4)-1 2*elements(j,4) ...
        2*elements(j,5)-1 2*elements(j,5) 2*elements(j,6)-1 2*elements(j,6) ...
        2*elements(j,7)-1 2*elements(j,7) 2*elements(j,8)-1 2*elements(j,8)];
    
    % We don't want temp <= 0, instead, put them to the end
    for tempi = 5:8
        if elements(j,tempi) == 0
            temp(2*tempi-1) = 2*(size(coordinates,1)+1)-1;
            temp(2*tempi)   = 2*(size(coordinates,1)+1);
        end
    end
    
    pt1 = eleNeighborIndexAndEdge(3,2);
    pt2 = eleNeighborIndexAndEdge(3,3);
    pointOfx = zeros(5,1); pointOfy = zeros(5,1);
    
    pointOfx(1) = gqpt1*0.5*(coordinates(pt2,1)-coordinates(pt1,1))+0.5*(coordinates(pt2,1)+coordinates(pt1,1));
    pointOfy(1) = gqpt1*0.5*(coordinates(pt2,2)-coordinates(pt1,2))+0.5*(coordinates(pt2,2)+coordinates(pt1,2));
    % pointOfx(2) = gqpt2*0.5*(coordinates(pt2,1)-coordinates(pt1,1))+0.5*(coordinates(pt2,1)+coordinates(pt1,1));
    pointOfy(2) = gqpt2*0.5*(coordinates(pt2,2)-coordinates(pt1,2))+0.5*(coordinates(pt2,2)+coordinates(pt1,2));
    % pointOfx(3) = gqpt3*0.5*(coordinates(pt2,1)-coordinates(pt1,1))+0.5*(coordinates(pt2,1)+coordinates(pt1,1));
    pointOfy(3) = gqpt3*0.5*(coordinates(pt2,2)-coordinates(pt1,2))+0.5*(coordinates(pt2,2)+coordinates(pt1,2));
    % pointOfx(4) = gqpt4*0.5*(coordinates(pt2,1)-coordinates(pt1,1))+0.5*(coordinates(pt2,1)+coordinates(pt1,1));
    pointOfy(4) = gqpt4*0.5*(coordinates(pt2,2)-coordinates(pt1,2))+0.5*(coordinates(pt2,2)+coordinates(pt1,2));
    % pointOfx(5) = gqpt5*0.5*(coordinates(pt2,1)-coordinates(pt1,1))+0.5*(coordinates(pt2,1)+coordinates(pt1,1));
    pointOfy(5) = gqpt5*0.5*(coordinates(pt2,2)-coordinates(pt1,2))+0.5*(coordinates(pt2,2)+coordinates(pt1,2));
    
    DUDX1 = zeros(size(pointOfx,1),size(pointOfx,1));
    DVDX1 = zeros(size(pointOfx,1),size(pointOfy,1));
    
    for tempi = 1 % :size(pointOfx,1)
        for tempj = 1: size(pointOfy,1)
            
            % ------ Calculate ksi and eta ------
            ksi = l(1)*pointOfx(tempi)*pointOfy(tempj) + l(2)*pointOfx(tempi) + l(3)*pointOfy(tempj) + l(4) ;
            eta = m(1)*pointOfx(tempi)*pointOfy(tempj) + m(2)*pointOfx(tempi) + m(3)*pointOfy(tempj) + m(4) ;
            
            % ------ Calculate N ------
            deltaPt5 = 0; deltaPt6 = 0; deltaPt7 = 0; deltaPt8 = 0;
            if elements(j,5) ~= 0
                deltaPt5 = 1;
            end
            if elements(j,6) ~= 0
                deltaPt6 = 1;
            end
            if elements(j,7) ~= 0
                deltaPt7 = 1;
            end
            if elements(j,8) ~= 0
                deltaPt8 = 1;
            end
            
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
            
            % ------- Calculate DU/Dx ---------
            tempDUDX = DN*U(temp);
            DUDX1(tempi,tempj) = tempDUDX(1);
            DVDX1(tempi,tempj) = tempDUDX(3);
            
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ========= Right side 1 element ===========
    
    elementtemp = eleNeighborIndexAndEdge(3,1);
    
    % ------ Find four corner points ------
    point1x = coordinates(elements(elementtemp,1),1);
    point1y = coordinates(elements(elementtemp,1),2);
    point2x = coordinates(elements(elementtemp,2),1);
    point2y = coordinates(elements(elementtemp,2),2);
    point3x = coordinates(elements(elementtemp,3),1);
    point3y = coordinates(elements(elementtemp,3),2);
    point4x = coordinates(elements(elementtemp,4),1);
    point4y = coordinates(elements(elementtemp,4),2);
    
    % ------ Find mid points 5/6/7/8 -------
    if elements(elementtemp,5) ~= 0
        point5x = coordinates(elements(elementtemp,5),1); 
        point5y = coordinates(elements(elementtemp,5),2);
    else
        point5x = 0; point5y = 0;
    end
    
    if elements(elementtemp,6) ~= 0
        point6x = coordinates(elements(elementtemp,6),1); 
        point6y = coordinates(elements(elementtemp,6),2);
    else
        point6x = 0; point6y = 0;
    end
    
    if elements(elementtemp,7) ~= 0
        point7x = coordinates(elements(elementtemp,7),1); 
        point7y = coordinates(elements(elementtemp,7),2);
    else
        point7x = 0; point7y = 0;
    end
    
    if elements(elementtemp,8) ~= 0
        point8x = coordinates(elements(elementtemp,8),1); 
        point8y = coordinates(elements(elementtemp,8),2);
    else
        point8x = 0; point8y = 0;
    end
    
    % ------ Calculate ksi and eta --------
    lMatrix = [point1x*point1y point1x point1y 1;
        point2x*point2y point2x point2y 1;
        point3x*point3y point3x point3y 1;
        point4x*point4y point4x point4y 1];
    
    lb = [-1;1;1;-1];
    l = linsolve(lMatrix,lb);
    
    mb = [-1;-1;1;1];
    m = linsolve(lMatrix,mb);
    
    % ------ Find the element nodal indices ------
    temp = [2*elements(elementtemp,1)-1 2*elements(elementtemp,1) ...
        2*elements(elementtemp,2)-1 2*elements(elementtemp,2) ...
        2*elements(elementtemp,3)-1 2*elements(elementtemp,3) ...
        2*elements(elementtemp,4)-1 2*elements(elementtemp,4) ...
        2*elements(elementtemp,5)-1 2*elements(elementtemp,5) ...
        2*elements(elementtemp,6)-1 2*elements(elementtemp,6) ...
        2*elements(elementtemp,7)-1 2*elements(elementtemp,7) ...
        2*elements(elementtemp,8)-1 2*elements(elementtemp,8)];
    
    % We don't want temp <= 0, instead, put them to the end
    for tempi = 5:8
        if elements(elementtemp,tempi) == 0
            temp(2*tempi-1) = 2*(size(coordinates,1)+1)-1;
            temp(2*tempi)   = 2*(size(coordinates,1)+1);
        end
    end
    
    
    DUDXRight1 = zeros(size(pointOfx,1),size(pointOfy,1));
    DVDXRight1 = zeros(size(pointOfx,1),size(pointOfy,1));
    
    for tempi = 1  %:size(pointOfx,1)
        for tempj = 1:size(pointOfy,1)
            
            % ------ Calculate ksi and eta -------
            ksi = l(1)*pointOfx(tempi)*pointOfy(tempj) + l(2)*pointOfx(tempi) + l(3)*pointOfy(tempj) + l(4);
            eta = m(1)*pointOfx(tempi)*pointOfy(tempj) + m(2)*pointOfx(tempi) + m(3)*pointOfy(tempj) + m(4);
            
            % ------ Calculate N ------
            deltaPt5 = 0; deltaPt6 = 0; deltaPt7 = 0; deltaPt8 = 0;
            if elements(elementtemp,5) ~= 0
                deltaPt5 = 1;
            end
            if elements(elementtemp,6) ~= 0
                deltaPt6 = 1;
            end
            if elements(elementtemp,7) ~= 0
                deltaPt7 = 1;
            end
            if elements(elementtemp,8) ~= 0
                deltaPt8 = 1;
            end
            
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
            
            % ------- Calculate DU/Dx ---------
            tempDUDX = DN*U(temp);
            DUDXRight1(tempi,tempj) = tempDUDX(1);
            DVDXRight1(tempi,tempj) = tempDUDX(3);
            
        end
    end
    
end


% ======== Right side 2 neighbor =========
if eleNeighborIndexAndEdge(4,1) ~= 0
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ========= Element itself ================
    
    % ------ Find four corner points ------
    point1x = coordinates(elements(j,1),1);
    point1y = coordinates(elements(j,1),2);
    point2x = coordinates(elements(j,2),1);
    point2y = coordinates(elements(j,2),2);
    point3x = coordinates(elements(j,3),1);
    point3y = coordinates(elements(j,3),2);
    point4x = coordinates(elements(j,4),1);
    point4y = coordinates(elements(j,4),2);
    
    % ------ Find mid points 5/6/7/8 -------
    if elements(j,5) ~= 0
        point5x = coordinates(elements(j,5),1); point5y = coordinates(elements(j,5),2);
    else
        point5x = 0; point5y = 0;
    end
    
    if elements(j,6) ~= 0
        point6x = coordinates(elements(j,6),1); point6y = coordinates(elements(j,6),2);
    else
        point6x = 0; point6y = 0;
    end
    
    if elements(j,7) ~= 0
        point7x = coordinates(elements(j,7),1); point7y = coordinates(elements(j,7),2);
    else
        point7x = 0; point7y = 0;
    end
    
    if elements(j,8) ~= 0
        point8x = coordinates(elements(j,8),1); point8y = coordinates(elements(j,8),2);
    else
        point8x = 0; point8y = 0;
    end
    
    % ------ Calculate ksi and eta --------
    lMatrix = [point1x*point1y point1x point1y 1;
        point2x*point2y point2x point2y 1;
        point3x*point3y point3x point3y 1;
        point4x*point4y point4x point4y 1];
    
    lb = [-1;1;1;-1];
    l = linsolve(lMatrix,lb);
    
    mb = [-1;-1;1;1];
    m = linsolve(lMatrix,mb);
     
    % ------ Find the element nodal indices ------
    temp = [2*elements(j,1)-1 2*elements(j,1) 2*elements(j,2)-1 2*elements(j,2) ...
        2*elements(j,3)-1 2*elements(j,3) 2*elements(j,4)-1 2*elements(j,4) ...
        2*elements(j,5)-1 2*elements(j,5) 2*elements(j,6)-1 2*elements(j,6) ...
        2*elements(j,7)-1 2*elements(j,7) 2*elements(j,8)-1 2*elements(j,8)];
    
    % We don't want temp <= 0, instead, put them to the end
    for tempi = 5:8
        if elements(j,tempi) == 0
            temp(2*tempi-1) = 2*(size(coordinates,1)+1)-1;
            temp(2*tempi)   = 2*(size(coordinates,1)+1);
        end
    end
    
    pt1 = eleNeighborIndexAndEdge(4,2);
    pt2 = eleNeighborIndexAndEdge(4,3);
    pointOfx = zeros(5,1); pointOfy = zeros(5,1);
    
    % Following pointOfx(1:5) and pointOfy(1:5) are used for Gaussian Quadrature
    
    pointOfx(1) = gqpt1*0.5*(coordinates(pt2,1)-coordinates(pt1,1))+0.5*(coordinates(pt2,1)+coordinates(pt1,1));
    pointOfy(1) = gqpt1*0.5*(coordinates(pt2,2)-coordinates(pt1,2))+0.5*(coordinates(pt2,2)+coordinates(pt1,2));
    % pointOfx(2) = gqpt2*0.5*(coordinates(pt2,1)-coordinates(pt1,1))+0.5*(coordinates(pt2,1)+coordinates(pt1,1));
    pointOfy(2) = gqpt2*0.5*(coordinates(pt2,2)-coordinates(pt1,2))+0.5*(coordinates(pt2,2)+coordinates(pt1,2));
    % pointOfx(3) = gqpt3*0.5*(coordinates(pt2,1)-coordinates(pt1,1))+0.5*(coordinates(pt2,1)+coordinates(pt1,1));
    pointOfy(3) = gqpt3*0.5*(coordinates(pt2,2)-coordinates(pt1,2))+0.5*(coordinates(pt2,2)+coordinates(pt1,2));
    % pointOfx(4) = gqpt4*0.5*(coordinates(pt2,1)-coordinates(pt1,1))+0.5*(coordinates(pt2,1)+coordinates(pt1,1));
    pointOfy(4) = gqpt4*0.5*(coordinates(pt2,2)-coordinates(pt1,2))+0.5*(coordinates(pt2,2)+coordinates(pt1,2));
    % pointOfx(5) = gqpt5*0.5*(coordinates(pt2,1)-coordinates(pt1,1))+0.5*(coordinates(pt2,1)+coordinates(pt1,1));
    pointOfy(5) = gqpt5*0.5*(coordinates(pt2,2)-coordinates(pt1,2))+0.5*(coordinates(pt2,2)+coordinates(pt1,2));
    
    DUDX2 = zeros(size(pointOfx,1),size(pointOfy,1));
    DVDX2 = zeros(size(pointOfx,1),size(pointOfy,1));
    
    for tempi = 1 % :size(pointOfx,1)
        for tempj = 1: size(pointOfy,1)
            
            % ------ Calculate ksi and eta ------
            ksi = l(1)*pointOfx(tempi)*pointOfy(tempj) + l(2)*pointOfx(tempi) + l(3)*pointOfy(tempj) + l(4) ;
            eta = m(1)*pointOfx(tempi)*pointOfy(tempj) + m(2)*pointOfx(tempi) + m(3)*pointOfy(tempj) + m(4) ;
            
            % ------ Calculate N ------
            deltaPt5 = 0; deltaPt6 = 0; deltaPt7 = 0; deltaPt8 = 0;
            if elements(j,5) ~= 0
                deltaPt5 = 1;
            end
            if elements(j,6) ~= 0
                deltaPt6 = 1;
            end
            if elements(j,7) ~= 0
                deltaPt7 = 1;
            end
            if elements(j,8) ~= 0
                deltaPt8 = 1;
            end
            
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
            
            % ------- Calculate DU/Dx ---------
            tempDUDX = DN*U(temp);
            DUDX2(tempi,tempj) = tempDUDX(1);
            DVDX2(tempi,tempj) = tempDUDX(3);
            
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ========= Right side element ===========
    
    elementtemp = eleNeighborIndexAndEdge(4,1);
    
    % ------ Find four corner points ------
    point1x = coordinates(elements(elementtemp,1),1);
    point1y = coordinates(elements(elementtemp,1),2);
    point2x = coordinates(elements(elementtemp,2),1);
    point2y = coordinates(elements(elementtemp,2),2);
    point3x = coordinates(elements(elementtemp,3),1);
    point3y = coordinates(elements(elementtemp,3),2);
    point4x = coordinates(elements(elementtemp,4),1);
    point4y = coordinates(elements(elementtemp,4),2);
    
    % ------ Find mid points 5/6/7/8 -------
    if elements(elementtemp,5) ~= 0
        point5x = coordinates(elements(elementtemp,5),1); 
        point5y = coordinates(elements(elementtemp,5),2);
    else
        point5x = 0; point5y = 0;
    end
    
    if elements(elementtemp,6) ~= 0
        point6x = coordinates(elements(elementtemp,6),1); 
        point6y = coordinates(elements(elementtemp,6),2);
    else
        point6x = 0; point6y = 0;
    end
    
    if elements(elementtemp,7) ~= 0
        point7x = coordinates(elements(elementtemp,7),1); 
        point7y = coordinates(elements(elementtemp,7),2);
    else
        point7x = 0; point7y = 0;
    end
    
    if elements(elementtemp,8) ~= 0
        point8x = coordinates(elements(elementtemp,8),1); 
        point8y = coordinates(elements(elementtemp,8),2);
    else
        point8x = 0; point8y = 0;
    end
    
    % ------ Calculate ksi and eta --------
    lMatrix = [ point1x*point1y point1x point1y 1;
                point2x*point2y point2x point2y 1;
                point3x*point3y point3x point3y 1;
                point4x*point4y point4x point4y 1 ];
        
    lb = [-1;1;1;-1];
    l = linsolve(lMatrix,lb);
    
    mb = [-1;-1;1;1];
    m = linsolve(lMatrix,mb);
    
    % ------ Find the element nodal indices ------
    temp = [2*elements(elementtemp,1)-1 2*elements(elementtemp,1) ...
        2*elements(elementtemp,2)-1 2*elements(elementtemp,2) ...
        2*elements(elementtemp,3)-1 2*elements(elementtemp,3) ...
        2*elements(elementtemp,4)-1 2*elements(elementtemp,4) ...
        2*elements(elementtemp,5)-1 2*elements(elementtemp,5) ...
        2*elements(elementtemp,6)-1 2*elements(elementtemp,6) ...
        2*elements(elementtemp,7)-1 2*elements(elementtemp,7) ...
        2*elements(elementtemp,8)-1 2*elements(elementtemp,8)];
    
    % We don't want temp <= 0, instead, put them to the end
    for tempi = 5:8
        if elements(elementtemp,tempi) == 0
            temp(2*tempi-1) = 2*(size(coordinates,1)+1)-1;
            temp(2*tempi)   = 2*(size(coordinates,1)+1);
        end
    end
    
    
    DUDXRight2 = zeros(size(pointOfx,1),size(pointOfy,1));
    DVDXRight2 = zeros(size(pointOfx,1),size(pointOfy,1));
    
    for tempi = 1  %:size(pointOfx,1)
        for tempj = 1:size(pointOfy,1)
            
            % ------ Calculate ksi and eta -------
            ksi = l(1)*pointOfx(tempi)*pointOfy(tempj) + l(2)*pointOfx(tempi) + l(3)*pointOfy(tempj) + l(4);
            eta = m(1)*pointOfx(tempi)*pointOfy(tempj) + m(2)*pointOfx(tempi) + m(3)*pointOfy(tempj) + m(4);
            
            % ------ Calculate N ------
            deltaPt5 = 0; deltaPt6 = 0; deltaPt7 = 0; deltaPt8 = 0;
            if elements(elementtemp,5) ~= 0
                deltaPt5 = 1;
            end
            if elements(elementtemp,6) ~= 0
                deltaPt6 = 1;
            end
            if elements(elementtemp,7) ~= 0
                deltaPt7 = 1;
            end
            if elements(elementtemp,8) ~= 0
                deltaPt8 = 1;
            end
            
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
            
            % ------- Calculate DU/Dx ---------
            tempDUDX = DN*U(temp);
            DUDXRight2(tempi,tempj) = tempDUDX(1);
            DVDXRight2(tempi,tempj) = tempDUDX(3);
            
        end
    end
    
end


% ============ Compute jump error ===============

if eleNeighborIndexAndEdge(3,1) ~= 0
    if eleNeighborIndexAndEdge(4,1) ~= 0
        for tempi = 1:length(gqwt)
            for tempj = 1:length(gqwt)
                ErrEstJump = ErrEstJump + gqwt(tempi)*gqwt(tempj)*lengthOfElement*0.5*((DUDX1(tempi,tempj)-DUDXRight1(tempi,tempj))^2 + ...
                            (DVDX1(tempi,tempj)-DVDXRight1(tempi,tempj))^2) + ...
                            gqwt(tempj)*gqwt(tempi)*lengthOfElement*0.5*((DUDX2(tempi,tempj)-DUDXRight2(tempi,tempj))^2 + ...
                            (DVDX2(tempi,tempj)-DVDXRight2(tempi,tempj))^2);
            end
        end
    else
        for tempi = 1:length(gqwt)
            for tempj = 1:length(gqwt)
                ErrEstJump = ErrEstJump + gqwt(tempi)*gqwt(tempj)*lengthOfElement*((DUDX1(tempi,tempj)-DUDXRight1(tempi,tempj))^2 + ...
                             (DVDX1(tempi,tempj)-DVDXRight1(tempi,tempj))^2);
            end
        end
    end
end
 





%% % ======== Top side =========

if eleNeighborIndexAndEdge(5,1) ~= 0
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ========= Element itself ================
    
    % ------ Find four corner points ------
    point1x = coordinates(elements(j,1),1);
    point1y = coordinates(elements(j,1),2);
    point2x = coordinates(elements(j,2),1);
    point2y = coordinates(elements(j,2),2);
    point3x = coordinates(elements(j,3),1);
    point3y = coordinates(elements(j,3),2);
    point4x = coordinates(elements(j,4),1);
    point4y = coordinates(elements(j,4),2);
    
    lengthOfElement = point3x-point1x;
    
    % ------ Find mid points 5/6/7/8 -------
    if elements(j,5) ~= 0
        point5x = coordinates(elements(j,5),1); point5y = coordinates(elements(j,5),2);
    else
        point5x = 0; point5y = 0;
    end
    
    if elements(j,6) ~= 0
        point6x = coordinates(elements(j,6),1); point6y = coordinates(elements(j,6),2);
    else
        point6x = 0; point6y = 0;
    end
    
    if elements(j,7) ~= 0
        point7x = coordinates(elements(j,7),1); point7y = coordinates(elements(j,7),2);
    else
        point7x = 0; point7y = 0;
    end
    
    if elements(j,8) ~= 0
        point8x = coordinates(elements(j,8),1); point8y = coordinates(elements(j,8),2);
    else
        point8x = 0; point8y = 0;
    end
    
    % ------ Calculate ksi and eta --------
    lMatrix = [point1x*point1y point1x point1y 1;
        point2x*point2y point2x point2y 1;
        point3x*point3y point3x point3y 1;
        point4x*point4y point4x point4y 1];
    
    lb = [-1;1;1;-1];
    l = linsolve(lMatrix,lb);
    
    mb = [-1;-1;1;1];
    m = linsolve(lMatrix,mb);
     
    % ------ Find the element nodal indices ------
    temp = [2*elements(j,1)-1 2*elements(j,1) 2*elements(j,2)-1 2*elements(j,2) ...
        2*elements(j,3)-1 2*elements(j,3) 2*elements(j,4)-1 2*elements(j,4) ...
        2*elements(j,5)-1 2*elements(j,5) 2*elements(j,6)-1 2*elements(j,6) ...
        2*elements(j,7)-1 2*elements(j,7) 2*elements(j,8)-1 2*elements(j,8)];
    
    % We don't want temp <= 0, instead, put them to the end
    for tempi = 5:8
        if elements(j,tempi) == 0
            temp(2*tempi-1) = 2*(size(coordinates,1)+1)-1;
            temp(2*tempi)   = 2*(size(coordinates,1)+1);
        end
    end
    
    pt1 = eleNeighborIndexAndEdge(5,2);
    pt2 = eleNeighborIndexAndEdge(5,3);
    pointOfx = zeros(5,1); pointOfy = zeros(5,1);
    
    pointOfx(1) = gqpt1*0.5*(coordinates(pt2,1)-coordinates(pt1,1))+0.5*(coordinates(pt2,1)+coordinates(pt1,1));
    pointOfy(1) = gqpt1*0.5*(coordinates(pt2,2)-coordinates(pt1,2))+0.5*(coordinates(pt2,2)+coordinates(pt1,2));
    pointOfx(2) = gqpt2*0.5*(coordinates(pt2,1)-coordinates(pt1,1))+0.5*(coordinates(pt2,1)+coordinates(pt1,1));
    % pointOfy(2) = gqpt2*0.5*(coordinates(pt2,2)-coordinates(pt1,2))+0.5*(coordinates(pt2,2)+coordinates(pt1,2));
    pointOfx(3) = gqpt3*0.5*(coordinates(pt2,1)-coordinates(pt1,1))+0.5*(coordinates(pt2,1)+coordinates(pt1,1));
    % pointOfy(3) = gqpt3*0.5*(coordinates(pt2,2)-coordinates(pt1,2))+0.5*(coordinates(pt2,2)+coordinates(pt1,2));
    pointOfx(4) = gqpt4*0.5*(coordinates(pt2,1)-coordinates(pt1,1))+0.5*(coordinates(pt2,1)+coordinates(pt1,1));
    % pointOfy(4) = gqpt4*0.5*(coordinates(pt2,2)-coordinates(pt1,2))+0.5*(coordinates(pt2,2)+coordinates(pt1,2));
    pointOfx(5) = gqpt5*0.5*(coordinates(pt2,1)-coordinates(pt1,1))+0.5*(coordinates(pt2,1)+coordinates(pt1,1));
    % pointOfy(5) = gqpt5*0.5*(coordinates(pt2,2)-coordinates(pt1,2))+0.5*(coordinates(pt2,2)+coordinates(pt1,2));
    
    DUDY1 = zeros(size(pointOfx,1),size(pointOfy,1));
    DVDY1 = zeros(size(pointOfx,1),size(pointOfy,1));
    
    for tempi = 1:size(pointOfx,1)
        for tempj = 1 % :size(pointOfy,1)
            
            % ------ Calculate ksi and eta ------
            ksi = l(1)*pointOfx(tempi)*pointOfy(tempj) + l(2)*pointOfx(tempi) + l(3)*pointOfy(tempj) + l(4) ;
            eta = m(1)*pointOfx(tempi)*pointOfy(tempj) + m(2)*pointOfx(tempi) + m(3)*pointOfy(tempj) + m(4) ;
            
            % ------ Calculate N ------
            deltaPt5 = 0; deltaPt6 = 0; deltaPt7 = 0; deltaPt8 = 0;
            if elements(j,5) ~= 0
                deltaPt5 = 1;
            end
            if elements(j,6) ~= 0
                deltaPt6 = 1;
            end
            if elements(j,7) ~= 0
                deltaPt7 = 1;
            end
            if elements(j,8) ~= 0
                deltaPt8 = 1;
            end
            
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
            
            % ------- Calculate DU/Dx ---------
            tempDUDX = DN*U(temp);
            DUDY1(tempi,tempj) = tempDUDX(2);
            DVDY1(tempi,tempj) = tempDUDX(4);
            
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ========= Left side 1 element ===========
    
    elementtemp = eleNeighborIndexAndEdge(5,1);
    
    % ------ Find four corner points ------
    point1x = coordinates(elements(elementtemp,1),1);
    point1y = coordinates(elements(elementtemp,1),2);
    point2x = coordinates(elements(elementtemp,2),1);
    point2y = coordinates(elements(elementtemp,2),2);
    point3x = coordinates(elements(elementtemp,3),1);
    point3y = coordinates(elements(elementtemp,3),2);
    point4x = coordinates(elements(elementtemp,4),1);
    point4y = coordinates(elements(elementtemp,4),2);
    
    % ------ Find mid points 5/6/7/8 -------
    if elements(elementtemp,5) ~= 0
        point5x = coordinates(elements(elementtemp,5),1); 
        point5y = coordinates(elements(elementtemp,5),2);
    else
        point5x = 0; point5y = 0;
    end
    
    if elements(elementtemp,6) ~= 0
        point6x = coordinates(elements(elementtemp,6),1); 
        point6y = coordinates(elements(elementtemp,6),2);
    else
        point6x = 0; point6y = 0;
    end
    
    if elements(elementtemp,7) ~= 0
        point7x = coordinates(elements(elementtemp,7),1); 
        point7y = coordinates(elements(elementtemp,7),2);
    else
        point7x = 0; point7y = 0;
    end
    
    if elements(elementtemp,8) ~= 0
        point8x = coordinates(elements(elementtemp,8),1); 
        point8y = coordinates(elements(elementtemp,8),2);
    else
        point8x = 0; point8y = 0;
    end
    
    % ------ Calculate ksi and eta --------
    lMatrix = [point1x*point1y point1x point1y 1;
        point2x*point2y point2x point2y 1;
        point3x*point3y point3x point3y 1;
        point4x*point4y point4x point4y 1];
    
    lb = [-1;1;1;-1];
    l = linsolve(lMatrix,lb);
    
    mb = [-1;-1;1;1];
    m = linsolve(lMatrix,mb);
    
    % ------ Find the element nodal indices ------
    temp = [2*elements(elementtemp,1)-1 2*elements(elementtemp,1) ...
        2*elements(elementtemp,2)-1 2*elements(elementtemp,2) ...
        2*elements(elementtemp,3)-1 2*elements(elementtemp,3) ...
        2*elements(elementtemp,4)-1 2*elements(elementtemp,4) ...
        2*elements(elementtemp,5)-1 2*elements(elementtemp,5) ...
        2*elements(elementtemp,6)-1 2*elements(elementtemp,6) ...
        2*elements(elementtemp,7)-1 2*elements(elementtemp,7) ...
        2*elements(elementtemp,8)-1 2*elements(elementtemp,8)];
    
    % We don't want temp <= 0, instead, put them to the end
    for tempi = 5:8
        if elements(elementtemp,tempi) == 0
            temp(2*tempi-1) = 2*(size(coordinates,1)+1)-1;
            temp(2*tempi)   = 2*(size(coordinates,1)+1);
        end
    end
    
    
    DUDYTop1 = zeros(size(pointOfx,1),size(pointOfy,1));
    DVDYTop1 = zeros(size(pointOfx,1),size(pointOfy,1));
    
    for tempi = 1:size(pointOfx,1)
        for tempj = 1 %:size(pointOfy,1)
            
            % ------ Calculate ksi and eta -------
            ksi = l(1)*pointOfx(tempi)*pointOfy(tempj) + l(2)*pointOfx(tempi) + l(3)*pointOfy(tempj) + l(4);
            eta = m(1)*pointOfx(tempi)*pointOfy(tempj) + m(2)*pointOfx(tempi) + m(3)*pointOfy(tempj) + m(4);
            
            % ------ Calculate N ------
            deltaPt5 = 0; deltaPt6 = 0; deltaPt7 = 0; deltaPt8 = 0;
            if elements(elementtemp,5) ~= 0
                deltaPt5 = 1;
            end
            if elements(elementtemp,6) ~= 0
                deltaPt6 = 1;
            end
            if elements(elementtemp,7) ~= 0
                deltaPt7 = 1;
            end
            if elements(elementtemp,8) ~= 0
                deltaPt8 = 1;
            end
            
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
            
            % ------- Calculate DU/Dx ---------
            tempDUDX = DN*U(temp);
            DUDYTop1(tempi,tempj) = tempDUDX(2);
            DVDYTop1(tempi,tempj) = tempDUDX(4);
            
        end
    end
    
end


% ======== Left side 2 neighbor =========
if eleNeighborIndexAndEdge(6,1) ~= 0
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ========= Element itself ================
    
    % ------ Find four corner points ------
    point1x = coordinates(elements(j,1),1);
    point1y = coordinates(elements(j,1),2);
    point2x = coordinates(elements(j,2),1);
    point2y = coordinates(elements(j,2),2);
    point3x = coordinates(elements(j,3),1);
    point3y = coordinates(elements(j,3),2);
    point4x = coordinates(elements(j,4),1);
    point4y = coordinates(elements(j,4),2);
    
    % ------ Find mid points 5/6/7/8 -------
    if elements(j,5) ~= 0
        point5x = coordinates(elements(j,5),1); point5y = coordinates(elements(j,5),2);
    else
        point5x = 0; point5y = 0;
    end
    
    if elements(j,6) ~= 0
        point6x = coordinates(elements(j,6),1); point6y = coordinates(elements(j,6),2);
    else
        point6x = 0; point6y = 0;
    end
    
    if elements(j,7) ~= 0
        point7x = coordinates(elements(j,7),1); point7y = coordinates(elements(j,7),2);
    else
        point7x = 0; point7y = 0;
    end
    
    if elements(j,8) ~= 0
        point8x = coordinates(elements(j,8),1); point8y = coordinates(elements(j,8),2);
    else
        point8x = 0; point8y = 0;
    end
    
    % ------ Calculate ksi and eta --------
    lMatrix = [point1x*point1y point1x point1y 1;
        point2x*point2y point2x point2y 1;
        point3x*point3y point3x point3y 1;
        point4x*point4y point4x point4y 1];
    
    lb = [-1;1;1;-1];
    l = linsolve(lMatrix,lb);
    
    mb = [-1;-1;1;1];
    m = linsolve(lMatrix,mb);
     
    % ------ Find the element nodal indices ------
    temp = [2*elements(j,1)-1 2*elements(j,1) 2*elements(j,2)-1 2*elements(j,2) ...
        2*elements(j,3)-1 2*elements(j,3) 2*elements(j,4)-1 2*elements(j,4) ...
        2*elements(j,5)-1 2*elements(j,5) 2*elements(j,6)-1 2*elements(j,6) ...
        2*elements(j,7)-1 2*elements(j,7) 2*elements(j,8)-1 2*elements(j,8)];
    
    % We don't want temp <= 0, instead, put them to the end
    for tempi = 5:8
        if elements(j,tempi) == 0
            temp(2*tempi-1) = 2*(size(coordinates,1)+1)-1;
            temp(2*tempi)   = 2*(size(coordinates,1)+1);
        end
    end
    
    pt1 = eleNeighborIndexAndEdge(6,2);
    pt2 = eleNeighborIndexAndEdge(6,3);
    pointOfx = zeros(5,1); pointOfy = zeros(5,1);
    
    % Following pointOfx(1:5) and pointOfy(1:5) are used for Gaussian Quadrature
    
    pointOfx(1) = gqpt1*0.5*(coordinates(pt2,1)-coordinates(pt1,1))+0.5*(coordinates(pt2,1)+coordinates(pt1,1));
    pointOfy(1) = gqpt1*0.5*(coordinates(pt2,2)-coordinates(pt1,2))+0.5*(coordinates(pt2,2)+coordinates(pt1,2));
    pointOfx(2) = gqpt2*0.5*(coordinates(pt2,1)-coordinates(pt1,1))+0.5*(coordinates(pt2,1)+coordinates(pt1,1));
    % pointOfy(2) = gqpt2*0.5*(coordinates(pt2,2)-coordinates(pt1,2))+0.5*(coordinates(pt2,2)+coordinates(pt1,2));
    pointOfx(3) = gqpt3*0.5*(coordinates(pt2,1)-coordinates(pt1,1))+0.5*(coordinates(pt2,1)+coordinates(pt1,1));
    % pointOfy(3) = gqpt3*0.5*(coordinates(pt2,2)-coordinates(pt1,2))+0.5*(coordinates(pt2,2)+coordinates(pt1,2));
    pointOfx(4) = gqpt4*0.5*(coordinates(pt2,1)-coordinates(pt1,1))+0.5*(coordinates(pt2,1)+coordinates(pt1,1));
    % pointOfy(4) = gqpt4*0.5*(coordinates(pt2,2)-coordinates(pt1,2))+0.5*(coordinates(pt2,2)+coordinates(pt1,2));
    pointOfx(5) = gqpt5*0.5*(coordinates(pt2,1)-coordinates(pt1,1))+0.5*(coordinates(pt2,1)+coordinates(pt1,1));
    % pointOfy(5) = gqpt5*0.5*(coordinates(pt2,2)-coordinates(pt1,2))+0.5*(coordinates(pt2,2)+coordinates(pt1,2));
    
    DUDY2 = zeros(size(pointOfx,1),size(pointOfy,1));
    DVDY2 = zeros(size(pointOfx,1),size(pointOfy,1));
    
    for tempi = 1:size(pointOfx,1)
        for tempj = 1  % : size(pointOfy,1)
            
            % ------ Calculate ksi and eta ------
            ksi = l(1)*pointOfx(tempi)*pointOfy(tempj) + l(2)*pointOfx(tempi) + l(3)*pointOfy(tempj) + l(4) ;
            eta = m(1)*pointOfx(tempi)*pointOfy(tempj) + m(2)*pointOfx(tempi) + m(3)*pointOfy(tempj) + m(4) ;
            
            % ------ Calculate N ------
            deltaPt5 = 0; deltaPt6 = 0; deltaPt7 = 0; deltaPt8 = 0;
            if elements(j,5) ~= 0
                deltaPt5 = 1;
            end
            if elements(j,6) ~= 0
                deltaPt6 = 1;
            end
            if elements(j,7) ~= 0
                deltaPt7 = 1;
            end
            if elements(j,8) ~= 0
                deltaPt8 = 1;
            end
            
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
            
            % ------- Calculate DU/Dx ---------
            tempDUDX = DN*U(temp);
            DUDY2(tempi,tempj) = tempDUDX(2);
            DVDY2(tempi,tempj) = tempDUDX(4);
            
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ========= Left side element ===========
    
    elementtemp = eleNeighborIndexAndEdge(6,1);
    
    % ------ Find four corner points ------
    point1x = coordinates(elements(elementtemp,1),1);
    point1y = coordinates(elements(elementtemp,1),2);
    point2x = coordinates(elements(elementtemp,2),1);
    point2y = coordinates(elements(elementtemp,2),2);
    point3x = coordinates(elements(elementtemp,3),1);
    point3y = coordinates(elements(elementtemp,3),2);
    point4x = coordinates(elements(elementtemp,4),1);
    point4y = coordinates(elements(elementtemp,4),2);
    
    % ------ Find mid points 5/6/7/8 -------
    if elements(elementtemp,5) ~= 0
        point5x = coordinates(elements(elementtemp,5),1); 
        point5y = coordinates(elements(elementtemp,5),2);
    else
        point5x = 0; point5y = 0;
    end
    
    if elements(elementtemp,6) ~= 0
        point6x = coordinates(elements(elementtemp,6),1); 
        point6y = coordinates(elements(elementtemp,6),2);
    else
        point6x = 0; point6y = 0;
    end
    
    if elements(elementtemp,7) ~= 0
        point7x = coordinates(elements(elementtemp,7),1); 
        point7y = coordinates(elements(elementtemp,7),2);
    else
        point7x = 0; point7y = 0;
    end
    
    if elements(elementtemp,8) ~= 0
        point8x = coordinates(elements(elementtemp,8),1); 
        point8y = coordinates(elements(elementtemp,8),2);
    else
        point8x = 0; point8y = 0;
    end
    
    % ------ Calculate ksi and eta --------
    lMatrix = [ point1x*point1y point1x point1y 1;
                point2x*point2y point2x point2y 1;
                point3x*point3y point3x point3y 1;
                point4x*point4y point4x point4y 1 ];
        
    lb = [-1;1;1;-1];
    l = linsolve(lMatrix,lb);
    
    mb = [-1;-1;1;1];
    m = linsolve(lMatrix,mb);
    
    % ------ Find the element nodal indices ------
    temp = [2*elements(elementtemp,1)-1 2*elements(elementtemp,1) ...
        2*elements(elementtemp,2)-1 2*elements(elementtemp,2) ...
        2*elements(elementtemp,3)-1 2*elements(elementtemp,3) ...
        2*elements(elementtemp,4)-1 2*elements(elementtemp,4) ...
        2*elements(elementtemp,5)-1 2*elements(elementtemp,5) ...
        2*elements(elementtemp,6)-1 2*elements(elementtemp,6) ...
        2*elements(elementtemp,7)-1 2*elements(elementtemp,7) ...
        2*elements(elementtemp,8)-1 2*elements(elementtemp,8)];
    
    % We don't want temp <= 0, instead, put them to the end
    for tempi = 5:8
        if elements(elementtemp,tempi) == 0
            temp(2*tempi-1) = 2*(size(coordinates,1)+1)-1;
            temp(2*tempi)   = 2*(size(coordinates,1)+1);
        end
    end
    
    
    DUDYTop2 = zeros(size(pointOfx,1),size(pointOfy,1));
    DVDYTop2 = zeros(size(pointOfx,1),size(pointOfy,1));
    
    for tempi = 1:size(pointOfx,1)
        for tempj = 1  %:size(pointOfy,1)
            
            % ------ Calculate ksi and eta -------
            ksi = l(1)*pointOfx(tempi)*pointOfy(tempj) + l(2)*pointOfx(tempi) + l(3)*pointOfy(tempj) + l(4);
            eta = m(1)*pointOfx(tempi)*pointOfy(tempj) + m(2)*pointOfx(tempi) + m(3)*pointOfy(tempj) + m(4);
            
            % ------ Calculate N ------
            deltaPt5 = 0; deltaPt6 = 0; deltaPt7 = 0; deltaPt8 = 0;
            if elements(elementtemp,5) ~= 0
                deltaPt5 = 1;
            end
            if elements(elementtemp,6) ~= 0
                deltaPt6 = 1;
            end
            if elements(elementtemp,7) ~= 0
                deltaPt7 = 1;
            end
            if elements(elementtemp,8) ~= 0
                deltaPt8 = 1;
            end
            
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
            
            % ------- Calculate DU/Dx ---------
            tempDUDX = DN*U(temp);
            DUDYTop2(tempi,tempj) = tempDUDX(2);
            DVDYTop2(tempi,tempj) = tempDUDX(4);
            
        end
    end
    
end


% ============ Compute jump error ===============

if eleNeighborIndexAndEdge(5,1) ~= 0
    if eleNeighborIndexAndEdge(6,1) ~= 0
        for tempi = 1:length(gqwt)
            for tempj = 1:length(gqwt)
                ErrEstJump = ErrEstJump + gqwt(tempi)*gqwt(tempj)*lengthOfElement*0.5*((DUDY1(tempi,tempj)-DUDYTop1(tempi,tempj))^2 ...
                             + (DVDY1(tempi,tempj)-DVDYTop1(tempi,tempj))^2) + ...
                             gqwt(tempi)*gqwt(tempj)*lengthOfElement*0.5*((DUDY2(tempi,tempj)-DUDYTop2(tempi,tempj))^2 ...
                             + (DVDY2(tempi,tempj)-DVDYTop2(tempi,tempj))^2);
            end
        end
    else
        for tempi = 1:length(gqwt)
            for tempj = 1:length(gqwt)
                ErrEstJump = ErrEstJump + gqwt(tempi)*gqwt(tempj)*lengthOfElement*((DUDY1(tempi,tempj)-DUDYTop1(tempi,tempj))^2 ...
                             + (DVDY1(tempi,tempj)-DVDYTop1(tempi,tempj))^2);
            end
        end
    end
end


 

%% % ======== Bottom side =========

if eleNeighborIndexAndEdge(7,1) ~= 0
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ========= Element itself ================
    
    % ------ Find four corner points ------
    point1x = coordinates(elements(j,1),1);
    point1y = coordinates(elements(j,1),2);
    point2x = coordinates(elements(j,2),1);
    point2y = coordinates(elements(j,2),2);
    point3x = coordinates(elements(j,3),1);
    point3y = coordinates(elements(j,3),2);
    point4x = coordinates(elements(j,4),1);
    point4y = coordinates(elements(j,4),2);
    
    lengthOfElement = point3x-point1x;
    
    % ------ Find mid points 5/6/7/8 -------
    if elements(j,5) ~= 0
        point5x = coordinates(elements(j,5),1); point5y = coordinates(elements(j,5),2);
    else
        point5x = 0; point5y = 0;
    end
    
    if elements(j,6) ~= 0
        point6x = coordinates(elements(j,6),1); point6y = coordinates(elements(j,6),2);
    else
        point6x = 0; point6y = 0;
    end
    
    if elements(j,7) ~= 0
        point7x = coordinates(elements(j,7),1); point7y = coordinates(elements(j,7),2);
    else
        point7x = 0; point7y = 0;
    end
    
    if elements(j,8) ~= 0
        point8x = coordinates(elements(j,8),1); point8y = coordinates(elements(j,8),2);
    else
        point8x = 0; point8y = 0;
    end
    
    % ------ Calculate ksi and eta --------
    lMatrix = [point1x*point1y point1x point1y 1;
        point2x*point2y point2x point2y 1;
        point3x*point3y point3x point3y 1;
        point4x*point4y point4x point4y 1];
    
    lb = [-1;1;1;-1];
    l = linsolve(lMatrix,lb);
    
    mb = [-1;-1;1;1];
    m = linsolve(lMatrix,mb);
     
    % ------ Find the element nodal indices ------
    temp = [2*elements(j,1)-1 2*elements(j,1) 2*elements(j,2)-1 2*elements(j,2) ...
        2*elements(j,3)-1 2*elements(j,3) 2*elements(j,4)-1 2*elements(j,4) ...
        2*elements(j,5)-1 2*elements(j,5) 2*elements(j,6)-1 2*elements(j,6) ...
        2*elements(j,7)-1 2*elements(j,7) 2*elements(j,8)-1 2*elements(j,8)];
    
    % We don't want temp <= 0, instead, put them to the end
    for tempi = 5:8
        if elements(j,tempi) == 0
            temp(2*tempi-1) = 2*(size(coordinates,1)+1)-1;
            temp(2*tempi)   = 2*(size(coordinates,1)+1);
        end
    end
    
    pt1 = eleNeighborIndexAndEdge(7,2);
    pt2 = eleNeighborIndexAndEdge(7,3);
    pointOfx = zeros(5,1); pointOfy = zeros(5,1);
    
    pointOfx(1) = gqpt1*0.5*(coordinates(pt2,1)-coordinates(pt1,1))+0.5*(coordinates(pt2,1)+coordinates(pt1,1));
    pointOfy(1) = gqpt1*0.5*(coordinates(pt2,2)-coordinates(pt1,2))+0.5*(coordinates(pt2,2)+coordinates(pt1,2));
    pointOfx(2) = gqpt2*0.5*(coordinates(pt2,1)-coordinates(pt1,1))+0.5*(coordinates(pt2,1)+coordinates(pt1,1));
    % pointOfy(2) = gqpt2*0.5*(coordinates(pt2,2)-coordinates(pt1,2))+0.5*(coordinates(pt2,2)+coordinates(pt1,2));
    pointOfx(3) = gqpt3*0.5*(coordinates(pt2,1)-coordinates(pt1,1))+0.5*(coordinates(pt2,1)+coordinates(pt1,1));
    % pointOfy(3) = gqpt3*0.5*(coordinates(pt2,2)-coordinates(pt1,2))+0.5*(coordinates(pt2,2)+coordinates(pt1,2));
    pointOfx(4) = gqpt4*0.5*(coordinates(pt2,1)-coordinates(pt1,1))+0.5*(coordinates(pt2,1)+coordinates(pt1,1));
    % pointOfy(4) = gqpt4*0.5*(coordinates(pt2,2)-coordinates(pt1,2))+0.5*(coordinates(pt2,2)+coordinates(pt1,2));
    pointOfx(5) = gqpt5*0.5*(coordinates(pt2,1)-coordinates(pt1,1))+0.5*(coordinates(pt2,1)+coordinates(pt1,1));
    % pointOfy(5) = gqpt5*0.5*(coordinates(pt2,2)-coordinates(pt1,2))+0.5*(coordinates(pt2,2)+coordinates(pt1,2));
    
    DUDY1 = zeros(size(pointOfx,1),size(pointOfy,1));
    DVDY1 = zeros(size(pointOfx,1),size(pointOfy,1));
    
    for tempi = 1:size(pointOfx,1)
        for tempj = 1 % :size(pointOfy,1)
            
            % ------ Calculate ksi and eta ------
            ksi = l(1)*pointOfx(tempi)*pointOfy(tempj) + l(2)*pointOfx(tempi) + l(3)*pointOfy(tempj) + l(4) ;
            eta = m(1)*pointOfx(tempi)*pointOfy(tempj) + m(2)*pointOfx(tempi) + m(3)*pointOfy(tempj) + m(4) ;
            
            % ------ Calculate N ------
            deltaPt5 = 0; deltaPt6 = 0; deltaPt7 = 0; deltaPt8 = 0;
            if elements(j,5) ~= 0
                deltaPt5 = 1;
            end
            if elements(j,6) ~= 0
                deltaPt6 = 1;
            end
            if elements(j,7) ~= 0
                deltaPt7 = 1;
            end
            if elements(j,8) ~= 0
                deltaPt8 = 1;
            end
            
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
            
            % ------- Calculate DU/Dx ---------
            tempDUDX = DN*U(temp);
            DUDY1(tempi,tempj) = tempDUDX(2);
            DVDY1(tempi,tempj) = tempDUDX(4);
            
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ========= Bottom side 1 element ===========
    
    elementtemp = eleNeighborIndexAndEdge(7,1);
    
    % ------ Find four corner points ------
    point1x = coordinates(elements(elementtemp,1),1);
    point1y = coordinates(elements(elementtemp,1),2);
    point2x = coordinates(elements(elementtemp,2),1);
    point2y = coordinates(elements(elementtemp,2),2);
    point3x = coordinates(elements(elementtemp,3),1);
    point3y = coordinates(elements(elementtemp,3),2);
    point4x = coordinates(elements(elementtemp,4),1);
    point4y = coordinates(elements(elementtemp,4),2);
    
    % ------ Find mid points 5/6/7/8 -------
    if elements(elementtemp,5) ~= 0
        point5x = coordinates(elements(elementtemp,5),1); 
        point5y = coordinates(elements(elementtemp,5),2);
    else
        point5x = 0; point5y = 0;
    end
    
    if elements(elementtemp,6) ~= 0
        point6x = coordinates(elements(elementtemp,6),1); 
        point6y = coordinates(elements(elementtemp,6),2);
    else
        point6x = 0; point6y = 0;
    end
    
    if elements(elementtemp,7) ~= 0
        point7x = coordinates(elements(elementtemp,7),1); 
        point7y = coordinates(elements(elementtemp,7),2);
    else
        point7x = 0; point7y = 0;
    end
    
    if elements(elementtemp,8) ~= 0
        point8x = coordinates(elements(elementtemp,8),1); 
        point8y = coordinates(elements(elementtemp,8),2);
    else
        point8x = 0; point8y = 0;
    end
    
    % ------ Calculate ksi and eta --------
    lMatrix = [point1x*point1y point1x point1y 1;
        point2x*point2y point2x point2y 1;
        point3x*point3y point3x point3y 1;
        point4x*point4y point4x point4y 1];
    
    lb = [-1;1;1;-1];
    l = linsolve(lMatrix,lb);
    
    mb = [-1;-1;1;1];
    m = linsolve(lMatrix,mb);
    
    % ------ Find the element nodal indices ------
    temp = [2*elements(elementtemp,1)-1 2*elements(elementtemp,1) ...
        2*elements(elementtemp,2)-1 2*elements(elementtemp,2) ...
        2*elements(elementtemp,3)-1 2*elements(elementtemp,3) ...
        2*elements(elementtemp,4)-1 2*elements(elementtemp,4) ...
        2*elements(elementtemp,5)-1 2*elements(elementtemp,5) ...
        2*elements(elementtemp,6)-1 2*elements(elementtemp,6) ...
        2*elements(elementtemp,7)-1 2*elements(elementtemp,7) ...
        2*elements(elementtemp,8)-1 2*elements(elementtemp,8)];
    
    % We don't want temp <= 0, instead, put them to the end
    for tempi = 5:8
        if elements(elementtemp,tempi) == 0
            temp(2*tempi-1) = 2*(size(coordinates,1)+1)-1;
            temp(2*tempi)   = 2*(size(coordinates,1)+1);
        end
    end
    
    
    DUDYBottom1 = zeros(size(pointOfx,1),size(pointOfy,1));
    DVDYBottom1 = zeros(size(pointOfx,1),size(pointOfy,1));
    
    for tempi = 1:size(pointOfx,1)
        for tempj = 1 %:size(pointOfy,1)
            
            % ------ Calculate ksi and eta -------
            ksi = l(1)*pointOfx(tempi)*pointOfy(tempj) + l(2)*pointOfx(tempi) + l(3)*pointOfy(tempj) + l(4);
            eta = m(1)*pointOfx(tempi)*pointOfy(tempj) + m(2)*pointOfx(tempi) + m(3)*pointOfy(tempj) + m(4);
            
            % ------ Calculate N ------
            deltaPt5 = 0; deltaPt6 = 0; deltaPt7 = 0; deltaPt8 = 0;
            if elements(elementtemp,5) ~= 0
                deltaPt5 = 1;
            end
            if elements(elementtemp,6) ~= 0
                deltaPt6 = 1;
            end
            if elements(elementtemp,7) ~= 0
                deltaPt7 = 1;
            end
            if elements(elementtemp,8) ~= 0
                deltaPt8 = 1;
            end
            
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
            
            % ------- Calculate DU/Dx ---------
            tempDUDX = DN*U(temp);
            DUDYBottom1(tempi,tempj) = tempDUDX(2);
            DVDYBottom1(tempi,tempj) = tempDUDX(4);
            
        end
    end
    
end


% ======== Bottom side 2 neighbor =========
if eleNeighborIndexAndEdge(8,1) ~= 0
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ========= Element itself ================
    
    % ------ Find four corner points ------
    point1x = coordinates(elements(j,1),1);
    point1y = coordinates(elements(j,1),2);
    point2x = coordinates(elements(j,2),1);
    point2y = coordinates(elements(j,2),2);
    point3x = coordinates(elements(j,3),1);
    point3y = coordinates(elements(j,3),2);
    point4x = coordinates(elements(j,4),1);
    point4y = coordinates(elements(j,4),2);
    
    % ------ Find mid points 5/6/7/8 -------
    if elements(j,5) ~= 0
        point5x = coordinates(elements(j,5),1); point5y = coordinates(elements(j,5),2);
    else
        point5x = 0; point5y = 0;
    end
    
    if elements(j,6) ~= 0
        point6x = coordinates(elements(j,6),1); point6y = coordinates(elements(j,6),2);
    else
        point6x = 0; point6y = 0;
    end
    
    if elements(j,7) ~= 0
        point7x = coordinates(elements(j,7),1); point7y = coordinates(elements(j,7),2);
    else
        point7x = 0; point7y = 0;
    end
    
    if elements(j,8) ~= 0
        point8x = coordinates(elements(j,8),1); point8y = coordinates(elements(j,8),2);
    else
        point8x = 0; point8y = 0;
    end
    
    % ------ Calculate ksi and eta --------
    lMatrix = [point1x*point1y point1x point1y 1;
        point2x*point2y point2x point2y 1;
        point3x*point3y point3x point3y 1;
        point4x*point4y point4x point4y 1];
    
    lb = [-1;1;1;-1];
    l = linsolve(lMatrix,lb);
    
    mb = [-1;-1;1;1];
    m = linsolve(lMatrix,mb);
     
    % ------ Find the element nodal indices ------
    temp = [2*elements(j,1)-1 2*elements(j,1) 2*elements(j,2)-1 2*elements(j,2) ...
        2*elements(j,3)-1 2*elements(j,3) 2*elements(j,4)-1 2*elements(j,4) ...
        2*elements(j,5)-1 2*elements(j,5) 2*elements(j,6)-1 2*elements(j,6) ...
        2*elements(j,7)-1 2*elements(j,7) 2*elements(j,8)-1 2*elements(j,8)];
    
    % We don't want temp <= 0, instead, put them to the end
    for tempi = 5:8
        if elements(j,tempi) == 0
            temp(2*tempi-1) = 2*(size(coordinates,1)+1)-1;
            temp(2*tempi)   = 2*(size(coordinates,1)+1);
        end
    end
    
    pt1 = eleNeighborIndexAndEdge(8,2);
    pt2 = eleNeighborIndexAndEdge(8,3);
    pointOfx = zeros(5,1); pointOfy = zeros(5,1);
    
    % Following pointOfx(1:5) and pointOfy(1:5) are used for Gaussian Quadrature
    
    pointOfx(1) = gqpt1*0.5*(coordinates(pt2,1)-coordinates(pt1,1))+0.5*(coordinates(pt2,1)+coordinates(pt1,1));
    pointOfy(1) = gqpt1*0.5*(coordinates(pt2,2)-coordinates(pt1,2))+0.5*(coordinates(pt2,2)+coordinates(pt1,2));
    pointOfx(2) = gqpt2*0.5*(coordinates(pt2,1)-coordinates(pt1,1))+0.5*(coordinates(pt2,1)+coordinates(pt1,1));
    % pointOfy(2) = gqpt2*0.5*(coordinates(pt2,2)-coordinates(pt1,2))+0.5*(coordinates(pt2,2)+coordinates(pt1,2));
    pointOfx(3) = gqpt3*0.5*(coordinates(pt2,1)-coordinates(pt1,1))+0.5*(coordinates(pt2,1)+coordinates(pt1,1));
    % pointOfy(3) = gqpt3*0.5*(coordinates(pt2,2)-coordinates(pt1,2))+0.5*(coordinates(pt2,2)+coordinates(pt1,2));
    pointOfx(4) = gqpt4*0.5*(coordinates(pt2,1)-coordinates(pt1,1))+0.5*(coordinates(pt2,1)+coordinates(pt1,1));
    % pointOfy(4) = gqpt4*0.5*(coordinates(pt2,2)-coordinates(pt1,2))+0.5*(coordinates(pt2,2)+coordinates(pt1,2));
    pointOfx(5) = gqpt5*0.5*(coordinates(pt2,1)-coordinates(pt1,1))+0.5*(coordinates(pt2,1)+coordinates(pt1,1));
    % pointOfy(5) = gqpt5*0.5*(coordinates(pt2,2)-coordinates(pt1,2))+0.5*(coordinates(pt2,2)+coordinates(pt1,2));
    
    DUDY2 = zeros(size(pointOfx,1),size(pointOfy,1));
    DVDY2 = zeros(size(pointOfx,1),size(pointOfy,1));
    
    for tempi = 1:size(pointOfx,1)
        for tempj = 1  % : size(pointOfy,1)
            
            % ------ Calculate ksi and eta ------
            ksi = l(1)*pointOfx(tempi)*pointOfy(tempj) + l(2)*pointOfx(tempi) + l(3)*pointOfy(tempj) + l(4) ;
            eta = m(1)*pointOfx(tempi)*pointOfy(tempj) + m(2)*pointOfx(tempi) + m(3)*pointOfy(tempj) + m(4) ;
            
            % ------ Calculate N ------
            deltaPt5 = 0; deltaPt6 = 0; deltaPt7 = 0; deltaPt8 = 0;
            if elements(j,5) ~= 0
                deltaPt5 = 1;
            end
            if elements(j,6) ~= 0
                deltaPt6 = 1;
            end
            if elements(j,7) ~= 0
                deltaPt7 = 1;
            end
            if elements(j,8) ~= 0
                deltaPt8 = 1;
            end
            
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
            
            % ------- Calculate DU/Dx ---------
            tempDUDX = DN*U(temp);
            DUDY2(tempi,tempj) = tempDUDX(2);
            DVDY2(tempi,tempj) = tempDUDX(4);
            
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ========= Bottom side element ===========
    
    elementtemp = eleNeighborIndexAndEdge(8,1);
    
    % ------ Find four corner points ------
    point1x = coordinates(elements(elementtemp,1),1);
    point1y = coordinates(elements(elementtemp,1),2);
    point2x = coordinates(elements(elementtemp,2),1);
    point2y = coordinates(elements(elementtemp,2),2);
    point3x = coordinates(elements(elementtemp,3),1);
    point3y = coordinates(elements(elementtemp,3),2);
    point4x = coordinates(elements(elementtemp,4),1);
    point4y = coordinates(elements(elementtemp,4),2);
    
    % ------ Find mid points 5/6/7/8 -------
    if elements(elementtemp,5) ~= 0
        point5x = coordinates(elements(elementtemp,5),1); 
        point5y = coordinates(elements(elementtemp,5),2);
    else
        point5x = 0; point5y = 0;
    end
    
    if elements(elementtemp,6) ~= 0
        point6x = coordinates(elements(elementtemp,6),1); 
        point6y = coordinates(elements(elementtemp,6),2);
    else
        point6x = 0; point6y = 0;
    end
    
    if elements(elementtemp,7) ~= 0
        point7x = coordinates(elements(elementtemp,7),1); 
        point7y = coordinates(elements(elementtemp,7),2);
    else
        point7x = 0; point7y = 0;
    end
    
    if elements(elementtemp,8) ~= 0
        point8x = coordinates(elements(elementtemp,8),1); 
        point8y = coordinates(elements(elementtemp,8),2);
    else
        point8x = 0; point8y = 0;
    end
    
    % ------ Calculate ksi and eta --------
    lMatrix = [ point1x*point1y point1x point1y 1;
                point2x*point2y point2x point2y 1;
                point3x*point3y point3x point3y 1;
                point4x*point4y point4x point4y 1 ];
        
    lb = [-1;1;1;-1];
    l = linsolve(lMatrix,lb);
    
    mb = [-1;-1;1;1];
    m = linsolve(lMatrix,mb);
    
    % ------ Find the element nodal indices ------
    temp = [2*elements(elementtemp,1)-1 2*elements(elementtemp,1) ...
        2*elements(elementtemp,2)-1 2*elements(elementtemp,2) ...
        2*elements(elementtemp,3)-1 2*elements(elementtemp,3) ...
        2*elements(elementtemp,4)-1 2*elements(elementtemp,4) ...
        2*elements(elementtemp,5)-1 2*elements(elementtemp,5) ...
        2*elements(elementtemp,6)-1 2*elements(elementtemp,6) ...
        2*elements(elementtemp,7)-1 2*elements(elementtemp,7) ...
        2*elements(elementtemp,8)-1 2*elements(elementtemp,8)];
    
    % We don't want temp <= 0, instead, put them to the end
    for tempi = 5:8
        if elements(elementtemp,tempi) == 0
            temp(2*tempi-1) = 2*(size(coordinates,1)+1)-1;
            temp(2*tempi)   = 2*(size(coordinates,1)+1);
        end
    end
    
    
    DUDYBottom2 = zeros(size(pointOfx,1),size(pointOfy,1));
    DVDYBottom2 = zeros(size(pointOfx,1),size(pointOfy,1));
    
    for tempi = 1:size(pointOfx,1)
        for tempj = 1  %:size(pointOfy,1)
            
            % ------ Calculate ksi and eta -------
            ksi = l(1)*pointOfx(tempi)*pointOfy(tempj) + l(2)*pointOfx(tempi) + l(3)*pointOfy(tempj) + l(4);
            eta = m(1)*pointOfx(tempi)*pointOfy(tempj) + m(2)*pointOfx(tempi) + m(3)*pointOfy(tempj) + m(4);
            
            % ------ Calculate N ------
            deltaPt5 = 0; deltaPt6 = 0; deltaPt7 = 0; deltaPt8 = 0;
            if elements(elementtemp,5) ~= 0
                deltaPt5 = 1;
            end
            if elements(elementtemp,6) ~= 0
                deltaPt6 = 1;
            end
            if elements(elementtemp,7) ~= 0
                deltaPt7 = 1;
            end
            if elements(elementtemp,8) ~= 0
                deltaPt8 = 1;
            end
            
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
            
            % ------- Calculate DU/Dx ---------
            tempDUDX = DN*U(temp);
            DUDYBottom2(tempi,tempj) = tempDUDX(2);
            DVDYBottom2(tempi,tempj) = tempDUDX(4);
            
        end
    end
    
end


% ============ Compute jump error ===============

if eleNeighborIndexAndEdge(7,1) ~= 0
    if eleNeighborIndexAndEdge(8,1) ~= 0
        for tempi = 1:length(gqwt)
            for tempj = 1:length(gqwt)
                ErrEstJump = ErrEstJump + gqwt(tempi)*gqwt(tempj)*lengthOfElement*0.5*((DUDY1(tempi,tempj)-DUDYBottom1(tempi,tempj))^2 ...
                             + (DVDY1(tempi,tempj)-DVDYBottom1(tempi,tempj))^2) + ...
                             gqwt(tempi)*gqwt(tempj)*lengthOfElement*0.5*((DUDY2(tempi,tempj)-DUDYBottom2(tempi,tempj))^2 ...
                             + (DVDY2(tempi,tempj)-DVDYBottom2(tempi,tempj))^2);
            end
        end
    else
        for tempi = 1:length(gqwt)
            for tempj = 1:length(gqwt)
                ErrEstJump = ErrEstJump + gqwt(tempi)*gqwt(tempj)*lengthOfElement*((DUDY1(tempi,tempj)-DUDYBottom1(tempi,tempj))^2 ...
                             + (DVDY1(tempi,tempj)-DVDYBottom1(tempi,tempj))^2);
            end
        end
    end
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



    




















