function pointInTriangleOrNot = funPointInTriangleCheck(point1x,point1y,point2x,point2y,point3x,point3y,pointOfx,pointOfy)

p1x = pointOfx; p1y = pointOfy; p2x = point1x; p2y = point1y; p3x = point2x; p3y = point2y;
b1 = sign( (p1x-p3x)*(p2y-p3y) - (p2x-p3x)*(p1y-p3y) );

p1x = pointOfx; p1y = pointOfy; p2x = point2x; p2y = point2y; p3x = point3x; p3y = point3y;
b2 = sign( (p1x-p3x)*(p2y-p3y) - (p2x-p3x)*(p1y-p3y) );

p1x = pointOfx; p1y = pointOfy; p2x = point3x; p2y = point3y; p3x = point1x; p3y = point1y;
b3 = sign( (p1x-p3x)*(p2y-p3y) - (p2x-p3x)*(p1y-p3y) );

if (b1==b2) && (b2==b3)
    pointInTriangleOrNot = 1; % point is inside of the triangle
else
    pointInTriangleOrNot = 0; % point is outside of the triangle
end



