function ErrEstInt = aPostErrEstIntSq(coordinates,elements,j,f,g,Df,U)

DfDx = Df.DfDx; DfDy = Df.DfDy; DfDxStartx = Df.DfAxis(1); DfDxStarty = Df.DfAxis(3);

ErrEstInt = 0; % Initialize ErrEstInt

% ------ Find four corner points ------
point1x = coordinates(elements(j,1),1);
point1y = coordinates(elements(j,1),2);
point2x = coordinates(elements(j,2),1);
point2y = coordinates(elements(j,2),2);
point3x = coordinates(elements(j,3),1);
point3y = coordinates(elements(j,3),2);
point4x = coordinates(elements(j,4),1);
point4y = coordinates(elements(j,4),2);

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

for pointOfx = point1x:point3x
    for pointOfy = point1y:point3y
        
        % ------ Calculate ksi and eta ------
        ksi = l(1)*pointOfx*pointOfy + l(2)*pointOfx + l(3)*pointOfy + l(4) ;
        eta = m(1)*pointOfx*pointOfy + m(2)*pointOfx + m(3)*pointOfy + m(4) ;
        
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
        
        % ------ Here approximate Dg(x+u)=Df(x) ------
        Df = [DfDx(pointOfx-DfDxStartx, pointOfy-DfDxStarty);
            DfDy(pointOfx-DfDxStartx, pointOfy-DfDxStarty)];
        
        temp1 = [pointOfx;pointOfy]+N*U;
        tempErrEstInt = sum(((f(pointOfx ,pointOfy) - fungInterpolation_g(temp1(1), temp1(2), g(floor(temp1(1))-1:floor(temp1(1))+2, floor(temp1(2))-1:floor(temp1(2))+2))) * Df).^2);
        
        ErrEstInt = ErrEstInt + tempErrEstInt;
        
    end
    
end

                    
                    
                    