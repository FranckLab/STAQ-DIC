%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function getRefinementPatch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [refinementPatch] = getRefinementPatch(refineT,g,elementsFEM,coordinatesFEM,winsize)

refinementPatch = refineT; % Initialize as itself
 
gRefineT = g(refineT);
eleEdgeLengthRefineT = winsize/(2^(gRefineT-1));

coordxLowerBound = min(coordinatesFEM(:,1));
coordxUpperBound = max(coordinatesFEM(:,1));
coordyLowerBound = min(coordinatesFEM(:,2));
coordyUpperBound = max(coordinatesFEM(:,2));

% =========================================================================
% Generate gImg 
gImg = sparse(coordxUpperBound ,coordyUpperBound );
for tempi = 1:size(elementsFEM,1)
    gImg(coordinatesFEM(elementsFEM(tempi,4),1),coordinatesFEM(elementsFEM(tempi,4),2)) = g(tempi);
end
  

% Find coordinates of element refineT left top corner
coordxLT = coordinatesFEM(elementsFEM(refineT,4),1);
coordyLT = coordinatesFEM(elementsFEM(refineT,4),2);

%% =========================================================================
% Find element refineT neighbors whose generation number is smaller than g(refineT)
% ------- # 1 -------
if coordxLT-2*eleEdgeLengthRefineT >= coordxLowerBound
    if gImg(coordxLT-2*eleEdgeLengthRefineT, coordyLT) == gRefineT-1
        refinementPatchLength = length(refinementPatch);
        refinementPatch(refinementPatchLength+1) = find(elementsFEM(:,5)==elementsFEM(refineT,1));
    end
end
% ------- # 2 -------
if coordxLT-2*eleEdgeLengthRefineT >= coordxLowerBound && coordyLT+eleEdgeLengthRefineT <= coordyUpperBound
    if gImg(coordxLT-2*eleEdgeLengthRefineT, coordyLT+eleEdgeLengthRefineT) == gRefineT-1
        refinementPatchLength = length(refinementPatch);
        refinementPatch(refinementPatchLength+1) = find(elementsFEM(:,5)==elementsFEM(refineT,4));
    end
end
% ------- # 3 -------
if coordxLT-eleEdgeLengthRefineT>= coordxLowerBound && coordyLT+2*eleEdgeLengthRefineT <= coordyUpperBound
    if gImg(coordxLT-eleEdgeLengthRefineT, coordyLT+2*eleEdgeLengthRefineT) == gRefineT-1
        refinementPatchLength = length(refinementPatch);
        refinementPatch(refinementPatchLength+1) = find(elementsFEM(:,8)==elementsFEM(refineT,4));
    end
end
% ------- # 4 -------
if coordyLT+2*eleEdgeLengthRefineT <= coordyUpperBound
    if gImg(coordxLT , coordyLT+2*eleEdgeLengthRefineT) == gRefineT-1
        refinementPatchLength = length(refinementPatch);
        refinementPatch(refinementPatchLength+1) = find(elementsFEM(:,8)==elementsFEM(refineT,3));
    end
end
% ------- # 5 -------
if coordxLT+eleEdgeLengthRefineT < coordxUpperBound && coordyLT+eleEdgeLengthRefineT <= coordyUpperBound
    if gImg(coordxLT+eleEdgeLengthRefineT, coordyLT+eleEdgeLengthRefineT) == gRefineT-1
        refinementPatchLength = length(refinementPatch);
        refinementPatch(refinementPatchLength+1) = find(elementsFEM(:,7)==elementsFEM(refineT,3));
    end
end
% ------- # 6 -------
if coordxLT+eleEdgeLengthRefineT < coordxUpperBound
    if gImg(coordxLT+eleEdgeLengthRefineT, coordyLT) == gRefineT-1
        refinementPatchLength = length(refinementPatch);
        refinementPatch(refinementPatchLength+1) = find(elementsFEM(:,7)==elementsFEM(refineT,2));
    end
end
% ------- # 7 -------
if coordyLT-eleEdgeLengthRefineT > coordyLowerBound
    if gImg(coordxLT, coordyLT-eleEdgeLengthRefineT) == gRefineT-1
        refinementPatchLength = length(refinementPatch);
        refinementPatch(refinementPatchLength+1) = find(elementsFEM(:,6)==elementsFEM(refineT,2));
    end
end
% ------- # 8 -------
if coordxLT-eleEdgeLengthRefineT >= coordxLowerBound && coordyLT-eleEdgeLengthRefineT > coordyLowerBound
    if gImg(coordxLT-eleEdgeLengthRefineT, coordyLT-eleEdgeLengthRefineT) == gRefineT-1
        refinementPatchLength = length(refinementPatch);
        refinementPatch(refinementPatchLength+1) = find(elementsFEM(:,6)==elementsFEM(refineT,1));
    end
end


% =========================================================================
% Find element refineT neighbors whose generation number is equal to g(refineT)
% % ------- # 1 -------
% if coordxLT-eleEdgeLengthRefineT >= coordxLowerBound
%     if gImg(coordxLT-eleEdgeLengthRefineT, coordyLT) == gRefineT
%         refinementPatchLength = length(refinementPatch);
%         refinementPatch(refinementPatchLength+1) = find(elementsFEM(:,3)==elementsFEM(refineT,4) );
%     end
% end
% % ------- # 2 -------
% if coordyLT+eleEdgeLengthRefineT <= coordyUpperBound
%     if gImg(coordxLT, coordyLT+eleEdgeLengthRefineT) == gRefineT
%         refinementPatchLength = length(refinementPatch);
%         refinementPatch(refinementPatchLength+1) = find(elementsFEM(:,1)==elementsFEM(refineT,4) );
%     end
% end
% % ------- # 3 -------
% if coordxLT+eleEdgeLengthRefineT < coordxUpperBound
%     if gImg(coordxLT+eleEdgeLengthRefineT, coordyLT) == gRefineT
%         refinementPatchLength = length(refinementPatch);
%         refinementPatch(refinementPatchLength+1) = find(elementsFEM(:,4)==elementsFEM(refineT,3)  );
%     end
% end
% % ------- # 4 -------
% if coordyLT-eleEdgeLengthRefineT > coordyLowerBound
%     if gImg(coordxLT, coordyLT-eleEdgeLengthRefineT) == gRefineT
%         refinementPatchLength = length(refinementPatch);
%         refinementPatch(refinementPatchLength+1) = find(elementsFEM(:,4)==elementsFEM(refineT,1) );
%     end
% end


