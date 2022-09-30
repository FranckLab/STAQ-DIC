function [markInside4,markOutside4] = markCircleInside(coordinates,elements4,C,R)
% markCircleInside: marks all the triangluar or quadrilateral elements with
%                   diamter > h inside the circle with radius R and center C

if ~isempty(elements4)
    
    ele4CenCoordx =  0.25*(  coordinates(elements4(:,1),1)  + coordinates(elements4(:,2),1) + ...
        coordinates(elements4(:,3),1) + coordinates(elements4(:,4),1)  );
    ele4CenCoordy =  0.25*(  coordinates(elements4(:,1),2)  + coordinates(elements4(:,2),2) + ...
        coordinates(elements4(:,3),2) + coordinates(elements4(:,4),2)  );
    
    ele4CenDist2Circ = sqrt((ele4CenCoordx-C(1)).^2 + (ele4CenCoordy-C(2)).^2);
    
    markInside4 = find(ele4CenDist2Circ < R+2);
    
    markOutside4 = setdiff([1:1:size(elements4,1)]', markInside4);
    
else
    markInside4 = [];
    markOutside4 = [];
    
end