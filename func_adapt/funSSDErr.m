function [ErrSSD] = funSSDErr(coordinates,elements,f,g,U,imageCutIndexOfx,imageCutIndexOfy)

ErrSSD = zeros(size(elements,1),1); % ErrPenalty = zeros(size(elements,1),1);
% h = waitbar(0,'Happy everyday!');

% parpool(4)
parfor j = 1: size(elements,1)
    
    % waitbar(j/size(elements,1));
    
    % ======= Step 1: find four corner points ======
    point1x = coordinates(elements(j,1),1);
    point1y = coordinates(elements(j,1),2);
    point2x = coordinates(elements(j,2),1);
    point2y = coordinates(elements(j,2),2);
    point3x = coordinates(elements(j,3),1);
    point3y = coordinates(elements(j,3),2);
    point4x = coordinates(elements(j,4),1);
    point4y = coordinates(elements(j,4),2);
    pointx = [point1x,point2x,point3x,point4x];
    pointy = [point1y,point2y,point3y,point4y];
    [cornersPositionOfx,~] = sort(pointx);
    [cornersPositionOfy,~] = sort(pointy);
    
    % ====== Comments =======
    % We can just check points within [cornersPositionOfx(1):cornersPositionOfx(4),
    %                                  cornersPositionOfy(1):cornersPositionOfy(4)]
    C = f(cornersPositionOfx(1):cornersPositionOfx(4), cornersPositionOfy(1):cornersPositionOfy(4));
    D = 0 * C;
    
    
    % ======= calculate ksi and eta coefficients =======
    lMatrix = [point1x*point1y point1x point1y 1;
        point2x*point2y point2x point2y 1;
        point3x*point3y point3x point3y 1;
        point4x*point4y point4x point4y 1];  % mMatrix = lMatrix;
    lb = [-1;1;1;-1];
    l = linsolve(lMatrix,lb);
    mb = [-1;-1;1;1];
    m = linsolve(lMatrix,mb);
    
    
    ci1 = 1; cj1 = 1; % index to count main loop
    for pointOfx = cornersPositionOfx(1):cornersPositionOfx(4)
       for pointOfy = cornersPositionOfy(1):cornersPositionOfy(4)
            
            % ====== calculate ksi and eta ======
            ksi = l(1)*pointOfx*pointOfy + l(2)*pointOfx + l(3)*pointOfy + l(4) ;
            eta = m(1)*pointOfx*pointOfy + m(2)*pointOfx + m(3)*pointOfy + m(4) ;
            
            % ====== calculate N matrix =======
            N1 = (1-ksi)*(1-eta)*0.25;
            N2 = (1+ksi)*(1-eta)*0.25;
            N3 = (1+ksi)*(1+eta)*0.25;
            N4 = (1-ksi)*(1+eta)*0.25;
            N = [N1 0 N2 0 N3 0 N4 0;
                0 N1 0 N2 0 N3 0 N4];
            
            temp = [2*elements(j,1)-1 2*elements(j,1) 2*elements(j,2)-1 2*elements(j,2)...
                2*elements(j,3)-1 2*elements(j,3) 2*elements(j,4)-1 2*elements(j,4)] ;
            temp1 = [pointOfx;pointOfy] + N*U(temp);
            
            tempg = g(floor(temp1(1)+imageCutIndexOfx)-1: floor(temp1(1)+imageCutIndexOfx)+2, ...
                floor(temp1(2)+imageCutIndexOfy)-1: floor(temp1(2)+imageCutIndexOfy)+2);
            D(cj1, ci1) =  fungInterpolation_g(temp1(1)+imageCutIndexOfx,temp1(2)+imageCutIndexOfy, tempg);
            ci1=ci1+1;
        end
        ci1=1;
        cj1=cj1+1;
    end
    
    %figure, surf(C); 
   % figure, surf(D); pause;
    
    % ErrSSD(j) = corr2(C,D);
    % Cmean = mean(C(:)); Cstd = std(C(:),1); Dmean = mean(D(:)); Dstd = std(D(:),1);
    % ErrSSD(j) = sum(sum(((C-Cmean)/Cstd - (D-Dmean)/Dstd).^2));
    ErrSSD(j) = sum(sum(((C )  - (D ) ).^2));
     
end