function ErrEstInt = aPostErrEstInt_ALSq(coordinates,elements,j,beta,mu,U,F,W,v,imgPyramidUnit,CrackOrNot)
 
ErrEstInt = 0; % Initialize ErrEstInt
U = [U;0;0]; v = [v;0;0]; F = [F;0;0;0;0]; W = [W;0;0;0;0];
FMinusW1 = F(1:2:end)-W(1:2:end); FMinusW2 = F(2:2:end)-W(2:2:end); UMinusv = U-v;

% ====== Gaussian quadrature parameter ======
if CrackOrNot == 1
    gqpt1 = -1; gqpt2 = 1; gqpt = [gqpt1,gqpt2];
    gqwt1 = 1;  gqwt2 = 1; gqwt = [gqwt1,gqwt2];
else
    gqpt1 = 0; gqpt2 = 0.538469; gqpt3 = -0.538469; gqpt4 = 0.90618; gqpt5 = -0.90618;
    gqwt1 = 0.568889; gqwt2 = 0.478629; gqwt3 = 0.478629; gqwt4 = 0.236927; gqwt5 = 0.236927;
    gqwt = [gqwt1,gqwt2,gqwt3,gqwt4,gqwt5]; gqpt = [gqpt1,gqpt2,gqpt3,gqpt4,gqpt5];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ------ Find four corner points ------
point1x = coordinates(elements(j,1),1);point1y = coordinates(elements(j,1),2);
point2x = coordinates(elements(j,2),1);point2y = coordinates(elements(j,2),2);
point3x = coordinates(elements(j,3),1);point3y = coordinates(elements(j,3),2);
point4x = coordinates(elements(j,4),1);point4y = coordinates(elements(j,4),2);

% ------ Find mid points 5/6/7/8 -------
try
    if elements(j,5) ~= 0
        point5x = coordinates(elements(j,5),1); point5y = coordinates(elements(j,5),2);
    else, point5x = 0; point5y = 0; end
catch
    point5x = 0; point5y = 0;
end
try
    if elements(j,6) ~= 0
        point6x = coordinates(elements(j,6),1); point6y = coordinates(elements(j,6),2);
    else, point6x = 0; point6y = 0; end
catch 
    point6x = 0; point6y = 0;
end
try
    if elements(j,7) ~= 0
        point7x = coordinates(elements(j,7),1); point7y = coordinates(elements(j,7),2);
    else, point7x = 0; point7y = 0; end
catch
    point7x = 0; point7y = 0;
end
try
    if elements(j,8) ~= 0
        point8x = coordinates(elements(j,8),1); point8y = coordinates(elements(j,8),2);
    else, point8x = 0; point8y = 0; end
catch
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

% ------ Set Gauss points ------
pt1 = elements(j,1); pt2 = elements(j,3);
pointOfx = zeros(length(gqwt),1); pointOfy = zeros(length(gqwt),1);
for tempi = 1:length(gqwt)
    pointOfx(tempi) = gqpt(tempi)*0.5*(coordinates(pt2,1)-coordinates(pt1,1))+0.5*(coordinates(pt2,1)+coordinates(pt1,1));
    pointOfy(tempi) = gqpt(tempi)*0.5*(coordinates(pt2,2)-coordinates(pt1,2))+0.5*(coordinates(pt2,2)+coordinates(pt1,2));
end

for tempi = 1:size(pointOfx,1)
    for tempj = 1:size(pointOfy,1)
        
        % ------ Calculate ksi and eta ------
        ksi = l(1)*pointOfx(tempi)*pointOfy(tempj) + l(2)*pointOfx(tempi) + l(3)*pointOfy(tempj) + l(4) ;
        eta = m(1)*pointOfx(tempi)*pointOfy(tempj) + m(2)*pointOfx(tempi) + m(3)*pointOfy(tempj) + m(4) ;

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
         tempDUDX = DN*FMinusW1(temp);
         DivFMinusW1 = tempDUDX(1)+tempDUDX(4);
         tempDUDX = DN*FMinusW2(temp);
         DivFMinusW2 = tempDUDX(1)+tempDUDX(4);
         
         tempErrEstInt = sum((  beta*DivFMinusW1 + mu*(U(temp(1:2:end))-UMinusv(temp(1:2:end)))  ).^2) + ...
                         sum((  beta*DivFMinusW2 + mu*(U(temp(2:2:end))-UMinusv(temp(2:2:end)))  ).^2);
         
         ErrEstInt = ErrEstInt + (point3x-point1x)*(point3y-point1y)*gqwt(tempi)*gqwt(tempj)*tempErrEstInt;
        
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

                    
                    
                    