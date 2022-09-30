function [longestEdgeMidPtCoord,longestEdgeNo,longestEdgeLength,eleCase] = funFindLongestEdge(coordinates)


% ------ Find three corner points ------
point1x = coordinates(1,1);
point1y = coordinates(1,2);
point2x = coordinates(2,1);
point2y = coordinates(2,2);
point3x = coordinates(3,1);
point3y = coordinates(3,2);


edgeLength1 = sqrt((point1x-point2x)^2 + (point1y-point2y)^2);
edgeLength2 = sqrt((point2x-point3x)^2 + (point2y-point3y)^2);
edgeLength3 = sqrt((point3x-point1x)^2 + (point3y-point1y)^2);

[longestEdgeLength,longestEdgeNo] = max([edgeLength1;edgeLength2;edgeLength3]);

if longestEdgeNo == 1
    longestEdgeMidPtCoord = [0.5*(point1x+point2x), 0.5*(point1y+point2y)];
    longestEdgeOppPtCoord = [point3x,point3y];
elseif longestEdgeNo == 2
    longestEdgeMidPtCoord = [0.5*(point2x+point3x), 0.5*(point2y+point3y)];
    longestEdgeOppPtCoord = [point1x,point1y];
elseif longestEdgeNo == 3
    longestEdgeMidPtCoord = [0.5*(point3x+point1x), 0.5*(point3y+point1y)];
    longestEdgeOppPtCoord = [point2x,point2y];
end

angletemp = atan2(-longestEdgeOppPtCoord(2)+longestEdgeMidPtCoord(2), -longestEdgeOppPtCoord(1)+longestEdgeMidPtCoord(1));

if abs(angletemp - pi/4) < 1e-3
    eleCase = 1;
elseif abs(angletemp - 2*pi/4) < 1e-3
    eleCase = 2;
elseif abs(angletemp - 3*pi/4) < 1e-3
    eleCase = 3;
elseif (abs(angletemp - 4*pi/4) < 1e-3) || (abs(angletemp + 4*pi/4) < 1e-3)
    eleCase = 4;
elseif abs(angletemp + 3*pi/4) < 1e-3
    eleCase = 5;
elseif abs(angletemp + 2*pi/4) < 1e-3
    eleCase = 6;
elseif abs(angletemp + 1*pi/4) < 1e-3
    eleCase = 7;
elseif abs(angletemp - 0) < 1e-3
    eleCase = 8;
else
    disp('Something wrong with this element!')
    eleCase = 0;
end

