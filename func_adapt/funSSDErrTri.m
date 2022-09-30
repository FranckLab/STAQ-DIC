function [ErrSSD] = funSSDErrTri(coordinates,elements,f,g,U,imageCutIndexOfx,imageCutIndexOfy)

ErrSSD = zeros(size(elements,1),1); % ErrPenalty = zeros(size(elements,1),1);
% h = waitbar(0,'Happy everyday!');

% parpool(4)
parfor j = 1: size(elements,1)
    
    % waitbar(j/size(elements,1));
    
    % ======= Step 1: find four corner points ======
    pt1x = coordinates(elements(j,1),1);
    pt1y = coordinates(elements(j,1),2);
    pt2x = coordinates(elements(j,2),1);
    pt2y = coordinates(elements(j,2),2);
    pt3x = coordinates(elements(j,3),1);
    pt3y = coordinates(elements(j,3),2);
     
    pointx = [pt1x,pt2x,pt3x ];
    pointy = [pt1y,pt2y,pt3y ];
    [cornersPositionOfx,~] = sort(pointx);
    [cornersPositionOfy,~] = sort(pointy);
    
    % ------ Compute triangle area --------
    TriArea = det([1 pt1x pt1y; 1 pt2x pt2y; 1 pt3x pt3y]);
    
    % ====== Comments =======
    % We can just check points within [cornersPositionOfx(1):cornersPositionOfx(4),
    %                                  cornersPositionOfy(1):cornersPositionOfy(4)]
    C = f(cornersPositionOfx(1):cornersPositionOfx(3), cornersPositionOfy(1):cornersPositionOfy(3));
    D = C;
    
    % ------ Find the element nodal indices ------
    tempIndexU = [2*elements(j,1)-1 2*elements(j,1) 2*elements(j,2)-1 2*elements(j,2) 2*elements(j,3)-1 2*elements(j,3) ];
        
    
    ci1 = 1; cj1 = 1; % index to count main loop
    for ptOfx = cornersPositionOfx(1):cornersPositionOfx(3)
       for ptOfy = cornersPositionOfy(1):cornersPositionOfy(3)
            
           
           % Judge pt is inside triangle or not
         ptInTriangleOrNot = funPointInTriangleCheck(pt1x,pt1y,pt2x,pt2y,pt3x,pt3y,ptOfx,ptOfy);
                
         if ptInTriangleOrNot == 1
                
             % ------ Calculate N ------
             N1 = det([1 ptOfx ptOfy; 1 pt2x pt2y; 1 pt3x pt3y])/TriArea;
             N2 = det([1 ptOfx ptOfy; 1 pt3x pt3y; 1 pt1x pt1y])/TriArea;
             N3 = det([1 ptOfx ptOfy; 1 pt1x pt1y; 1 pt2x pt2y])/TriArea;
             
             N = [N1 0 N2 0 N3 0;
                 0 N1 0 N2 0 N3];
              
            temp1 = [ptOfx;ptOfy] + N*U(tempIndexU);
            
            tempg = g(floor(temp1(1)+imageCutIndexOfx)-1: floor(temp1(1)+imageCutIndexOfx)+2, ...
                floor(temp1(2)+imageCutIndexOfy)-1: floor(temp1(2)+imageCutIndexOfy)+2);
            D(cj1, ci1) =  fungInterpolation_g(temp1(1)+imageCutIndexOfx,temp1(2)+imageCutIndexOfy, tempg);
            
         end
         
         ci1=ci1+1;
        end
        ci1=1;
        cj1=cj1+1;
    end
    
    % ErrSSD(j) = corr2(C,D);
    % Cmean = mean(C(:)); Cstd = std(C(:),1); Dmean = mean(D(:)); Dstd = std(D(:),1);
    % ErrSSD(j) = sum(sum(((C-Cmean)/Cstd - (D-Dmean)/Dstd).^2));
    ErrSSD(j) = sum(sum(((C )  - (D ) ).^2));
     
end