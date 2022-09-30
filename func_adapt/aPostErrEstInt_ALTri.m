function ErrEstInt = aPostErrEstInt_ALTri(coordinates,elements,j,beta,mu,U,F,W,v,imgPyramidUnit,CrackOrNot)

ErrEstInt = 0; % Initialize ErrEstInt
U = [U;0;0]; v = [v;0;0]; F = [F;0;0;0;0]; W = [W;0;0;0;0];
FMinusW1 = F(1:2:end)-W(1:2:end); FMinusW2 = F(2:2:end)-W(2:2:end); UMinusv = U-v;
 
% ------ Find three corner points ------
point1x = coordinates(elements(j,1),1);point1y = coordinates(elements(j,1),2);
point2x = coordinates(elements(j,2),1);point2y = coordinates(elements(j,2),2);
point3x = coordinates(elements(j,3),1);point3y = coordinates(elements(j,3),2);
 
% ------ Compute triangle area --------
TriArea =  det([1 point1x point1y; 1 point2x point2y; 1 point3x point3y]);
 
% ------ Calculate DN Matrix for CST ------
funDN1x = 1/(2*TriArea)*(point2y-point3y); funDN1y = 1/(2*TriArea)*(point3x-point2x);
funDN2x = 1/(2*TriArea)*(point3y-point1y); funDN2y = 1/(2*TriArea)*(point1x-point3x);
funDN3x = 1/(2*TriArea)*(point1y-point2y); funDN3y = 1/(2*TriArea)*(point2x-point1x);

DN = [ funDN1x 0 funDN2x 0 funDN3x 0  ;
       funDN1y 0 funDN2y 0 funDN3y 0  ;
       0 funDN1x 0 funDN2x 0 funDN3x  ;
       0 funDN1y 0 funDN2y 0 funDN3y  ];
         
% ------ Find the element nodal indices ------
temp = [2*elements(j,1)-1 2*elements(j,1) 2*elements(j,2)-1 2*elements(j,2) ...
        2*elements(j,3)-1 2*elements(j,3)];
 
for pointOfx = 1/3*(point1x+point2x+point3x)
    for pointOfy = 1/3*(point1y+point2y+point3y)
        
        % Judge point is inside triangle or not
        pointInTriangleOrNot = funPointInTriangleCheck(point1x,point1y,point2x,point2y,point3x,point3y,pointOfx,pointOfy);
        
        if pointInTriangleOrNot == 1
            
            % ------ Calculate N ------
            N1 = det([1 pointOfx pointOfy; 1 point2x point2y; 1 point3x point3y])/TriArea;
            N2 = det([1 pointOfx pointOfy; 1 point3x point3y; 1 point1x point1y])/TriArea;
            N3 = det([1 pointOfx pointOfy; 1 point1x point1y; 1 point2x point2y])/TriArea;
            
            N = [N1 0 N2 0 N3 0;
                0 N1 0 N2 0 N3];
            
            % ------- Calculate DU/Dx ---------
            tempDUDX = DN*FMinusW1(temp);
            DivFMinusW1 = tempDUDX(1)+tempDUDX(4);
            tempDUDX = DN*FMinusW2(temp);
            DivFMinusW2 = tempDUDX(1)+tempDUDX(4);
            
            tempErrEstInt = sum((  beta*DivFMinusW1 + mu*(U(temp(1:2:end))-UMinusv(temp(1:2:end)))  ).^2) + ...
                sum((  beta*DivFMinusW2 + mu*(U(temp(2:2:end))-UMinusv(temp(2:2:end)))  ).^2);
            
            ErrEstInt = ErrEstInt + TriArea*tempErrEstInt;
            
        end
    end
    
end
