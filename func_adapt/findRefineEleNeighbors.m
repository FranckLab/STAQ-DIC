function eleNeighbor = findRefineEleNeighbors(elementsFEM,elej,whichSide)
 
switch whichSide
    
    % ============================= case 5 ================================
    case 5
        
        % ----- Find right side element neighbor ------
        [row1,~] = find(elementsFEM == elementsFEM(elej,2));
        [row2,~] = find(elementsFEM == elementsFEM(elej,3));
        
        row3 = intersect(row1,row2);
        
        if length(row3) == 1
            % it's at the boundary, no element neighbor
             eleNeighbor = 0;
        else
            for tempi = 1:length(row3)
                if (row3(tempi)==elej)
                    continue;
                else
                    eleNeighbor = row3(tempi);
                end
            end
        end
        
    % ============================= case 6 ===============================
    case 6
        
        % ----- Find right side element neighbor ------
        [row1,~] = find(elementsFEM == elementsFEM(elej,3));
        [row2,~] = find(elementsFEM == elementsFEM(elej,4));
        
        row3 = intersect(row1,row2);
        
        if length(row3) == 1
            % it's at the boundary, no element neighbor
             eleNeighbor = 0;
        else
            for tempi = 1:length(row3)
                if (row3(tempi)==elej)
                    continue;
                else
                    eleNeighbor = row3(tempi);
                end
            end
        end
        
    % ============================= case 7 ===============================
    case 7
        
        % ----- Find right side element neighbor ------
        [row1,~] = find(elementsFEM == elementsFEM(elej,1));
        [row2,~] = find(elementsFEM == elementsFEM(elej,4));
        
        row3 = intersect(row1,row2);
        
        if length(row3) == 1
            % it's at the boundary, no element neighbor
             eleNeighbor = 0;
        else
            for tempi = 1:length(row3)
                if (row3(tempi)==elej)
                    continue;
                else
                    eleNeighbor = row3(tempi);
                end
            end
        end
        
    % ============================= case 8 ===============================
    case 8
        
        % ----- Find right side element neighbor ------
        [row1,~] = find(elementsFEM == elementsFEM(elej,1));
        [row2,~] = find(elementsFEM == elementsFEM(elej,2));
        
        row3 = intersect(row1,row2);
        
        if length(row3) == 1
            % it's at the boundary, no element neighbor
             eleNeighbor = 0;
        else
            for tempi = 1:length(row3)
                if (row3(tempi)==elej)
                    continue;
                else
                    eleNeighbor = row3(tempi);
                end
            end
        end
        
end





