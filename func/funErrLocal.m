function [SSDErr,f2] = funErrLocal(coordinates,elements,coordinatesFEM,f,g,U,F)

SSDErr = zeros(size(elements,1),1); f2 = 0*f;
h = waitbar(0);
for j = 1:size(elements ,1)
    waitbar(j/size(elements,1));
    x0 = coordinatesFEM(j,1); y0 = coordinatesFEM(j,2);
    
    point1 = elements(j,1); point3 = elements(j,3); 
    x1 = coordinates(point1,1); y1 = coordinates(point1,2);
    x3 = coordinates(point3,1); y3 = coordinates(point3,2);
    
    % Find region for f
    tempf = f(ceil(x1):floor(x3), ceil(y1):floor(y3));
    
    % Find region for g
    tempg = zeros(size(tempf,1), size(tempf,2));
     
    P = [F(4*j-3:4*j);U(2*j-1:2*j)];
    for tempi = ceil(x1) :floor(x3) 
        for tempj = ceil(y1) :floor(y3) 
            Wp = [x0;y0]+[1+P(1)  P(3) P(5) ; P(2) 1+P(4)  P(6) ]*[tempi-x0;tempj-y0;1];
            u22 = Wp(1); v22 = Wp(2);
            tempg(tempi+1-(ceil(x1)),tempj+1-(ceil(y1)))=...
                fungInterpolation_g(u22, v22, g(floor(u22)-1:floor(u22)+2, floor(v22)-1:floor(v22)+2));
        end
    end
    f2(ceil(x1) :floor(x3) , ceil(y1) :floor(y3) ) = tempg;
    
    SSDErr(j) = sum((tempf(:)-tempg(:)).^2);

end