function ErrEstInt = aPostErrEstInt(coordinates,f,g,Df,U)

DfDx = Df.DfDx; DfDy = Df.DfDy; DfDxStartx = Df.DfAxis(1); DfDxStarty = Df.DfAxis(3);

ErrEstInt = 0;

% ------ Find three corner pts ------
pt1x = coordinates(1,1);
pt1y = coordinates(1,2);
pt2x = coordinates(2,1);
pt2y = coordinates(2,2);
pt3x = coordinates(3,1);
pt3y = coordinates(3,2);

% ------ Compute triangle area --------
TriArea = det([1 pt1x pt1y; 1 pt2x pt2y; 1 pt3x pt3y]);

 

for ptOfx = min([pt1x,pt2x,pt3x]):max([pt1x,pt2x,pt3x])
    for ptOfy = min([pt1y,pt2y,pt3y]):max([pt1y,pt2y,pt3y])
        
        % Judge pt is inside triangle or not
        ptInTriangleOrNot = funPointInTriangleCheck(pt1x,pt1y,pt2x,pt2y,pt3x,pt3y,ptOfx,ptOfy);
        
        if ptInTriangleOrNot == 1
            
            % ------ Calculate N ------
            N1 = det([1 ptOfx ptOfy; 1 pt2x pt2y; 1 pt3x pt3y])/TriArea;
            N2 = det([1 ptOfx ptOfy; 1 pt3x pt3y; 1 pt1x pt1y])/TriArea;
            N3 = det([1 ptOfx ptOfy; 1 pt1x pt1y; 1 pt2x pt2y])/TriArea;
            
            N = [N1 0 N2 0 N3 0;
                0 N1 0 N2 0 N3];
            
            % ------ Here approximate Dg(x+u)=Df(x) ------
            Df = [DfDx(ptOfx-DfDxStartx, ptOfy-DfDxStarty);
                DfDy(ptOfx-DfDxStartx, ptOfy-DfDxStarty)];
            
            temp1 = [ptOfx;ptOfy]+N*U;
            tempErrEstInt = sum(((f(ptOfx ,ptOfy) - fungInterpolation_g(temp1(1), temp1(2), g(floor(temp1(1))-1:floor(temp1(1))+2, floor(temp1(2))-1:floor(temp1(2))+2))) * Df).^2);
            % tempErrEstInt = sum(((f(ptOfx ,ptOfy) - fungInterpolation_g(temp1(1), temp1(2), g(floor(temp1(1))-1:floor(temp1(1))+2, floor(temp1(2))-1:floor(temp1(2))+2))) ).^2);
            ErrEstInt = ErrEstInt + tempErrEstInt;
            
        end
        
    end
end
                    
                    
                    