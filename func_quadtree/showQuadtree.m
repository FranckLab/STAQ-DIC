function h2=showQuadtree(elements4,coordinates,u)
 

minCoordStep = 3;

xList = [ceil(min(coordinates(:,1))) : minCoordStep : floor(max(coordinates(:,1)))]';
yList = [ceil(min(coordinates(:,2))) : minCoordStep : floor(max(coordinates(:,2)))]';

% uGrid = gridfit( coordinates(:,1), coordinates(:,2), u, xList, yList ,'regularizer','springs'); uGrid=uGrid';
[xGrid,yGrid] = ndgrid(xList, yList);
F_u = scatteredInterpolant( coordinates(:,1), coordinates(:,2), u );
uGrid = F_u(xGrid,yGrid);


% Recover the image mask
ImgRefMask = zeros( xList(end), yList(end));
for eleInd = 1:size(elements4,1)
    xCoordMin = ceil(min(coordinates( elements4(eleInd,:),1 )));
    xCoordMax = floor(max(coordinates( elements4(eleInd,:),1 )));
    yCoordMin = ceil(min(coordinates( elements4(eleInd,:),2 )));
    yCoordMax = floor(max(coordinates( elements4(eleInd,:),2 )));
    ImgRefMask( xCoordMin:xCoordMax,  yCoordMin:yCoordMax ) = 1;
end
% testpts = [xGrid(:),yGrid(:)];
% xyz = [ coordinates( elements4(eleInd,1:4),1 ),  coordinates( elements4(eleInd,1:4), 2 )   ];
% in = inhull(testpts,xyz);

ImgRefMask(ImgRefMask==0) = nan;

% Compute the uGridmask
uGridMask = ImgRefMask(xList, yList);




h2=surf( xGrid, yGrid, uGrid.*uGridMask, 'EdgeColor','none','LineStyle','none' );

