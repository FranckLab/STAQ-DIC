function eleNeighborIndexAndEdge = findEleNeighbors(elementsFEM,coordinatesFEM,elementT,eleGeneration,winsize)

eleNeighborIndexAndEdge = zeros(8,3);

coordxLowerBound = min(coordinatesFEM(:,1));
coordxUpperBound = max(coordinatesFEM(:,1));
coordyLowerBound = min(coordinatesFEM(:,2));
coordyUpperBound = max(coordinatesFEM(:,2));
 
% Generate gImg = eleGeneration Image
gImg = sparse(coordxUpperBound , coordyUpperBound );
for tempi = 1:size(elementsFEM,1)
    gImg(coordinatesFEM(elementsFEM(tempi,4),1),coordinatesFEM(elementsFEM(tempi,4),2)) = eleGeneration(tempi);
end

% Find eleGeneration number of "j" element
gElementT = eleGeneration(elementT);
eleEdgeLengthElementT = winsize/(2^(gElementT-1));

% Find coordinates of element refineT left top corner
coordxLT = coordinatesFEM(elementsFEM(elementT,4),1);
coordyLT = coordinatesFEM(elementsFEM(elementT,4),2);


%% =========================================================================
% Find element elementT neighbors whose generation number is larger than eleGeneration(elementT)
if gElementT>1 
    
    % ------- # 1 -------
    if coordxLT-2*eleEdgeLengthElementT >= coordxLowerBound
        if gImg(coordxLT-2*eleEdgeLengthElementT, coordyLT) == gElementT-1
            elementtemp = find(elementsFEM(:,5)==elementsFEM(elementT,1));
            eleNeighborIndexAndEdge(1,:) = [ elementtemp, elementsFEM(elementtemp,3), elementsFEM(elementtemp,5)];
        end
    end
    % ------- # 2 -------
    if coordxLT-2*eleEdgeLengthElementT >= coordxLowerBound && coordyLT+eleEdgeLengthElementT <= coordyUpperBound
        if gImg(coordxLT-2*eleEdgeLengthElementT, coordyLT+eleEdgeLengthElementT) == gElementT-1
            elementtemp = find(elementsFEM(:,5)==elementsFEM(elementT,4));
            eleNeighborIndexAndEdge(1,:) = [ elementtemp, elementsFEM(elementtemp,5), elementsFEM(elementtemp,2)];
        end
    end
    % ------- # 3 -------
    if coordxLT-eleEdgeLengthElementT>= coordxLowerBound && coordyLT+2*eleEdgeLengthElementT <= coordyUpperBound
        if gImg(coordxLT-eleEdgeLengthElementT, coordyLT+2*eleEdgeLengthElementT) == gElementT-1
            elementtemp = find(elementsFEM(:,8)==elementsFEM(elementT,4));
            eleNeighborIndexAndEdge(3,:) = [ elementtemp, elementsFEM(elementtemp,2), elementsFEM(elementtemp,8)];
        end
    end
    % ------- # 4 -------
    if coordyLT+2*eleEdgeLengthElementT <= coordyUpperBound
        if gImg(coordxLT , coordyLT+2*eleEdgeLengthElementT) == gElementT-1
            elementtemp = find(elementsFEM(:,8)==elementsFEM(elementT,3));
            eleNeighborIndexAndEdge(3,:) = [ elementtemp, elementsFEM(elementtemp,8), elementsFEM(elementtemp,1)];
        end
    end
    % ------- # 5 -------
    if coordxLT+eleEdgeLengthElementT < coordxUpperBound && coordyLT+eleEdgeLengthElementT <= coordyUpperBound
        if gImg(coordxLT+eleEdgeLengthElementT, coordyLT+eleEdgeLengthElementT) == gElementT-1
            elementtemp = find(elementsFEM(:,7)==elementsFEM(elementT,3));
            eleNeighborIndexAndEdge(5,:) = [ elementtemp, elementsFEM(elementtemp,1), elementsFEM(elementtemp,7)];
        end
    end
    % ------- # 6 -------
    if coordxLT+eleEdgeLengthElementT < coordxUpperBound
        if gImg(coordxLT+eleEdgeLengthElementT, coordyLT) == gElementT-1
            elementtemp = find(elementsFEM(:,7)==elementsFEM(elementT,2));
            eleNeighborIndexAndEdge(5,:) = [ elementtemp, elementsFEM(elementtemp,7), elementsFEM(elementtemp,4)];
        end
    end
    % ------- # 7 -------
    if coordyLT-eleEdgeLengthElementT > coordyLowerBound
        if gImg(coordxLT, coordyLT-eleEdgeLengthElementT) == gElementT-1
            elementtemp = find(elementsFEM(:,6)==elementsFEM(elementT,2));
            eleNeighborIndexAndEdge(7,:) = [ elementtemp, elementsFEM(elementtemp,4), elementsFEM(elementtemp,6)];
        end
    end
    % ------- # 8 -------
    if coordxLT-eleEdgeLengthElementT >= coordxLowerBound && coordyLT-eleEdgeLengthElementT > coordyLowerBound
        if gImg(coordxLT-eleEdgeLengthElementT, coordyLT-eleEdgeLengthElementT) == gElementT-1
            elementtemp = find(elementsFEM(:,6)== elementsFEM(elementT,1));
            eleNeighborIndexAndEdge(7,:) = [ elementtemp, elementsFEM(elementtemp,6), elementsFEM(elementtemp,3)];
        end
    end
    
end

%% =========================================================================
% Find element refineT neighbors whose generation number is equal to g(refineT)
% ------- # 1 -------
if coordxLT-eleEdgeLengthElementT >= coordxLowerBound
    if gImg(coordxLT-eleEdgeLengthElementT, coordyLT) == gElementT
        elementtemp = find(elementsFEM(:,3)==elementsFEM(elementT,4));
        eleNeighborIndexAndEdge(1,:) = [ elementtemp, elementsFEM(elementtemp,3), elementsFEM(elementtemp,2)];
    end
end
% ------- # 2 -------
if coordyLT+eleEdgeLengthElementT <= coordyUpperBound
    if gImg(coordxLT, coordyLT+eleEdgeLengthElementT) == gElementT
        elementtemp = find(elementsFEM(:,1)==elementsFEM(elementT,4));
        eleNeighborIndexAndEdge(3,:) = [ elementtemp, elementsFEM(elementtemp,2), elementsFEM(elementtemp,1)];
    end
end
% ------- # 3 -------
if coordxLT+eleEdgeLengthElementT < coordxUpperBound
    if gImg(coordxLT+eleEdgeLengthElementT, coordyLT) == gElementT
        elementtemp = find(elementsFEM(:,4)==elementsFEM(elementT,3));
        eleNeighborIndexAndEdge(5,:) = [ elementtemp, elementsFEM(elementtemp,1), elementsFEM(elementtemp,4)];
    end
end
% ------- # 4 -------
if coordyLT-eleEdgeLengthElementT > coordyLowerBound
    if gImg(coordxLT, coordyLT-eleEdgeLengthElementT) == gElementT
        elementtemp = find(elementsFEM(:,4)==elementsFEM(elementT,1));
        eleNeighborIndexAndEdge(7,:) = [ elementtemp, elementsFEM(elementtemp,4), elementsFEM(elementtemp,3)];
    end
end


%% =========================================================================
% Find element refineT neighbors whose generation number is smaller to g(refineT)

% ------- # 1 & 2 -------
if coordxLT-0.5*eleEdgeLengthElementT >= coordxLowerBound 
    if coordyLT-0.5*eleEdgeLengthElementT > coordyLowerBound
        if gImg(coordxLT-0.5*eleEdgeLengthElementT, coordyLT-0.5*eleEdgeLengthElementT) == gElementT+1
            elementtemp = find(elementsFEM(:,2)==elementsFEM(elementT,1));
            eleNeighborIndexAndEdge(1,:) = [ elementtemp, elementsFEM(elementtemp,3), elementsFEM(elementtemp,2)];
        end
    end
    if gImg(coordxLT-0.5*eleEdgeLengthElementT,coordyLT) == gElementT+1
        elementtemp = find(elementsFEM(:,3)==elementsFEM(elementT,4));
        eleNeighborIndexAndEdge(2,:) = [ elementtemp, elementsFEM(elementtemp,3), elementsFEM(elementtemp,2)];
    end
end
 
% ------- # 3 & 4 -------
if coordyLT+0.5*eleEdgeLengthElementT <= coordyUpperBound
    if gImg(coordxLT, coordyLT+0.5*eleEdgeLengthElementT) == gElementT+1
        elementtemp = find(elementsFEM(:,1)==elementsFEM(elementT,4));
        eleNeighborIndexAndEdge(3,:) = [ elementtemp, elementsFEM(elementtemp,2), elementsFEM(elementtemp,1)];
    end
    if coordxLT+0.5*eleEdgeLengthElementT < coordxUpperBound
        if gImg(coordxLT+0.5*eleEdgeLengthElementT, coordyLT+0.5*eleEdgeLengthElementT) == gElementT+1
            elementtemp = find(elementsFEM(:,2)==elementsFEM(elementT,3));
            eleNeighborIndexAndEdge(4,:) = [ elementtemp, elementsFEM(elementtemp,2), elementsFEM(elementtemp,1)];
        end
    end
end

% ------- # 5 & 6 -------
if coordxLT+eleEdgeLengthElementT < coordxUpperBound 
    if gImg(coordxLT+eleEdgeLengthElementT, coordyLT) == gElementT+1
        elementtemp = find(elementsFEM(:,4)==elementsFEM(elementT,3));
        eleNeighborIndexAndEdge(5,:) = [ elementtemp, elementsFEM(elementtemp,1), elementsFEM(elementtemp,4)];
    end
    if coordyLT-0.5*eleEdgeLengthElementT > coordyLowerBound
        if gImg(coordxLT+eleEdgeLengthElementT, coordyLT-0.5*eleEdgeLengthElementT) == gElementT+1
            elementtemp = find(elementsFEM(:,1)==elementsFEM(elementT,2));
            eleNeighborIndexAndEdge(6,:) = [ elementtemp, elementsFEM(elementtemp,1), elementsFEM(elementtemp,4)];
        end
    end
end

% ------- # 7 & 8 -------
if coordyLT-eleEdgeLengthElementT > coordyLowerBound
    if coordxLT+0.5*eleEdgeLengthElementT < coordxUpperBound
        if gImg(coordxLT+0.5*eleEdgeLengthElementT, coordyLT-eleEdgeLengthElementT) == gElementT+1
            elementtemp = find(elementsFEM(:,3)==elementsFEM(elementT,2));
            eleNeighborIndexAndEdge(7,:) = [ elementtemp, elementsFEM(elementtemp,4), elementsFEM(elementtemp,3)];
        end
    end
    if gImg(coordxLT,coordyLT-eleEdgeLengthElementT) == gElementT+1
        elementtemp = find(elementsFEM(:,4)==elementsFEM(elementT,1));
        eleNeighborIndexAndEdge(8,:) = [ elementtemp, elementsFEM(elementtemp,4), elementsFEM(elementtemp,3)];
    end
end
 



