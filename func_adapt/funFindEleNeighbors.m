function [eleNeighbor2,eleNeighbor3,eleNeighbor4] = funFindEleNeighbors(elementsFEM,j,longestEdgeNo)
 
% ------ Find eleNeighbor4 ------
[row1,~] = find(elementsFEM == elementsFEM(j,1+mod(longestEdgeNo+2,3)));
[row2,~] = find(elementsFEM == elementsFEM(j,1+mod(longestEdgeNo+3,3)));

row3 = intersect(row1,row2);
if length(row3)==1
    % it's at the boundary, no sharing edge elements
    eleNeighbor4 = 0;
else
    for tempi = 1:2
        if (row3(tempi)==j)
            continue
        else
            eleNeighbor4 = row3(tempi);
        end
    end
end

 

% ------ Find eleNeighbor2 ------
[row1,~] = find(elementsFEM == elementsFEM(j,1+mod(longestEdgeNo+3,3)));
[row2,~] = find(elementsFEM == elementsFEM(j,1+mod(longestEdgeNo+1,3)));

row3 = intersect(row1,row2);
if length(row3)==1
    % it's at the boundary, no sharing edge elements
    eleNeighbor2 = 0;
else
    for tempi = 1:2
        if (row3(tempi)==j)
            continue
        else
            eleNeighbor2 = row3(tempi);
        end
    end
end

% ------ Find eleNeighbor3 ------
[row1,~] = find(elementsFEM == elementsFEM(j,1+mod(longestEdgeNo+1,3)));
[row2,~] = find(elementsFEM == elementsFEM(j,1+mod(longestEdgeNo+2,3)));

row3 = intersect(row1,row2);
if length(row3)==1
    % it's at the boundary, no sharing edge elements
    eleNeighbor3 = 0;
else
    for tempi = 1:2
        if (row3(tempi)==j)
            continue
        else
            eleNeighbor3 = row3(tempi);
        end
    end
end
    
