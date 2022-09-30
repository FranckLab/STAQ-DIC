function [coordinatesFEMNew,elementsFEMNew,UNew,FNew,refinedEleIDList,dirichletNew,neumannNew] = ...
    refineRecursiveTri(coordinatesFEM, elementsFEM, U,F, refinementPatch, refinementPatchAll, FinishRefineRecursive)

% Comment: for each single element refinement

% disp('==Just into refinerecursive function ===')
% FinishRefineRecursive
% refinementPatch
% disp('======')

coordinatesFEMNew = coordinatesFEM; elementsFEMNew = elementsFEM;
UNew = U; FNew = F;  dirichletNew = []; neumannNew = [];

while FinishRefineRecursive ~= 0
    
    refinedEleIDList = 0;
    
    FinishRefineRecursive = 0;
    refinementPatchtemp = zeros(2*length(refinementPatch),1);
    gRefinePatchtemp = zeros(2*length(refinementPatch),1);
    
    % ======= Get refinement patch R(J,T) =======
    % for tempi = 1:length(refinementPatch)
    for tempi = 1:1
        
        if refinementPatch(tempi) ~= 0
             
            % Find element longest edge length
            [~,longestEdgeNo,longestEdgeLength,~] = funFindLongestEdge(coordinatesFEMNew(elementsFEMNew(refinementPatch(tempi),:),:));
             
            % Find neighboring elements connected with the longest edge
            [~,~,eleNeighbor4] = funFindEleNeighbors(elementsFEMNew,refinementPatch(tempi),longestEdgeNo);
            
            if eleNeighbor4 ~= 0 % Not at the boundary
                
                [~,~,longestEdgeLengthNeighbor4,~] = funFindLongestEdge(coordinatesFEMNew(elementsFEMNew(eleNeighbor4,:),:));
                
                if (longestEdgeLength < longestEdgeLengthNeighbor4 ) 
                     
                    refinementPatchtemp = eleNeighbor4;
                     
                    FinishRefineRecursive = FinishRefineRecursive + 1;
                    
                    [coordinatesFEMNew,elementsFEMNew,UNew,FNew,refinedEleIDList,dirichletNew,neumannNew] = refineRecursiveTri(coordinatesFEM, elementsFEM, U,F, refinementPatchtemp, refinementPatchAll, FinishRefineRecursive) ;
                    
                    
                    % disp('=== FinishRefineRecursive is ===')
                    % FinishRefineRecursive
                    % FinishRefineRecursive = 0;
                    
                else
                    
                    FinishRefineRecursive = 0;
                     
                end
                
            else
            % At the boundary: elements can be directly refined
            % Do nothing, and it will automatically go to the following
            % bisection part.
                
               
                
            end
            
            
        end
        
    end
    
end

% disp('======== INTO triangulation part =========')

% ======= Triangulation mesh element bisection part =======
if FinishRefineRecursive == 0
 
    coordinatesFEM = coordinatesFEMNew; elementsFEM = elementsFEMNew; U= UNew;
    
    % ======= Get refinement patch R(J,T) =======
    for tempi = 1:1
        
          
        if refinementPatch(tempi) ~= 0
            
            % Find element longest edge length
            [longestEdgeMidPtCoord,longestEdgeNo,~,~] = funFindLongestEdge(coordinatesFEM(elementsFEM(refinementPatch(tempi),:),:));
            
             
            % Find neighboring elements connected with the longest edge
            [~,~,eleNeighbor4] = funFindEleNeighbors(elementsFEM,refinementPatch(tempi),longestEdgeNo);
             
            
            refinementPatchtemp(2*tempi-1) = refinementPatch(tempi);
            refinementPatchtemp(2*tempi) = eleNeighbor4;
            
            
            % Add one more mid-point after the end of coordinatesFEM
            LengthOfCoordsOld = size(coordinatesFEM,1);
            coordinatesFEMNew(LengthOfCoordsOld+1,:) = longestEdgeMidPtCoord;
             
            
            % Cover original elementsFEM by a new bisection element + add another bisection elements
            LengthOfElementsOld = size(elementsFEM,1);
            
 
            if eleNeighbor4 ~= 0 % Longest edge is not at the boundary, bisect both the element and the eleNeighbor4
                 
                % But we still need to know whether eleNeighbor4 belongs to
                % refinementPatchAll or not
                % checkeleNgb4InRefinePatchOrNot = ismember(eleNeighbor4,refinementPatchAll);
                
                % Find eleNeighbor4's longest edge No
                [~,eleNg4longestEdgeNo,~,~] = funFindLongestEdge(coordinatesFEM(elementsFEM(eleNeighbor4,:),:));
                
                  
                switch eleNg4longestEdgeNo
                    
                    case 1
                        elementsFEMNew(LengthOfElementsOld+2,:) = [elementsFEM(refinementPatchtemp(2*tempi),3), elementsFEM(refinementPatchtemp(2*tempi),1), LengthOfCoordsOld+1];
                        % if checkeleNgb4InRefinePatchOrNot == 0
                        %     elementsFEMNew(LengthOfElementsOld+3,:) = [elementsFEM(refinementPatchtemp(2*tempi),2:3), LengthOfCoordsOld+1];
                        % else
                            elementsFEMNew(eleNeighbor4,:) = [elementsFEM(refinementPatchtemp(2*tempi),2:3), LengthOfCoordsOld+1];
                        % end
                        % Interpolate displacements
                        UNew(2*LengthOfCoordsOld+1) = 0.5*(U(2*(elementsFEM(refinementPatchtemp(2*tempi),1))-1)+ ...
                                            U(2*(elementsFEM(refinementPatchtemp(2*tempi),2))-1));
                        UNew(2*LengthOfCoordsOld+2) = 0.5*(U(2*(elementsFEM(refinementPatchtemp(2*tempi),1)))+ ...
                                                U(2*(elementsFEM(refinementPatchtemp(2*tempi),2))));
                        FNew(4*LengthOfCoordsOld+1) = 0.5*(F(4*(elementsFEM(refinementPatchtemp(2*tempi),1))-3)+ ...
                                            F(4*(elementsFEM(refinementPatchtemp(2*tempi),2))-3));
                        FNew(4*LengthOfCoordsOld+2) = 0.5*(F(4*(elementsFEM(refinementPatchtemp(2*tempi),1))-2)+ ...
                                           F(4*(elementsFEM(refinementPatchtemp(2*tempi),2))-2));
                        FNew(4*LengthOfCoordsOld+3) = 0.5*(F(4*(elementsFEM(refinementPatchtemp(2*tempi),1))-1)+ ...
                                            F(4*(elementsFEM(refinementPatchtemp(2*tempi),2))-1));
                        FNew(4*LengthOfCoordsOld+4) = 0.5*(F(4*(elementsFEM(refinementPatchtemp(2*tempi),1)))+ ...
                                                F(4*(elementsFEM(refinementPatchtemp(2*tempi),2))));
                            
        
                    case 2
                        elementsFEMNew(LengthOfElementsOld+2,:) = [elementsFEM(refinementPatchtemp(2*tempi),3), elementsFEM(refinementPatchtemp(2*tempi),1), LengthOfCoordsOld+1];
                        % if checkeleNgb4InRefinePatchOrNot == 0
                        %     elementsFEMNew(LengthOfElementsOld+3,:) = [elementsFEM(refinementPatchtemp(2*tempi),1:2), LengthOfCoordsOld+1];
                        % else
                            elementsFEMNew(eleNeighbor4,:) = [elementsFEM(refinementPatchtemp(2*tempi),1:2), LengthOfCoordsOld+1];
                        % end
                        % Interpolate displacements
                        UNew(2*LengthOfCoordsOld+1) = 0.5*(U(2*(elementsFEM(refinementPatchtemp(2*tempi),2))-1)+U(2*(elementsFEM(refinementPatchtemp(2*tempi),3))-1));
                        UNew(2*LengthOfCoordsOld+2) = 0.5*(U(2*(elementsFEM(refinementPatchtemp(2*tempi),2)))+U(2*(elementsFEM(refinementPatchtemp(2*tempi),3))));                        
                        
                        FNew(4*LengthOfCoordsOld+1) = 0.5*(F(4*(elementsFEM(refinementPatchtemp(2*tempi),2))-3)+ ...
                                            F(4*(elementsFEM(refinementPatchtemp(2*tempi),3))-3));
                        FNew(4*LengthOfCoordsOld+2) = 0.5*(F(4*(elementsFEM(refinementPatchtemp(2*tempi),2))-2)+ ...
                                            F(4*(elementsFEM(refinementPatchtemp(2*tempi),3))-2));
                        FNew(4*LengthOfCoordsOld+3) = 0.5*(F(4*(elementsFEM(refinementPatchtemp(2*tempi),2))-1)+ ...
                                            F(4*(elementsFEM(refinementPatchtemp(2*tempi),3))-1));
                        FNew(4*LengthOfCoordsOld+4) = 0.5*(F(4*(elementsFEM(refinementPatchtemp(2*tempi),2)))+ ...
                                               F(4*(elementsFEM(refinementPatchtemp(2*tempi),3))));

                    case 3
                        elementsFEMNew(LengthOfElementsOld+2,:) = [elementsFEM(refinementPatchtemp(2*tempi),1:2), LengthOfCoordsOld+1];
                        % if checkeleNgb4InRefinePatchOrNot == 0
                        %     elementsFEMNew(LengthOfElementsOld+3,:) = [elementsFEM(refinementPatchtemp(2*tempi),2:3), LengthOfCoordsOld+1];
                        % else
                            elementsFEMNew(eleNeighbor4,:) = [elementsFEM(refinementPatchtemp(2*tempi),2:3), LengthOfCoordsOld+1];
                        % end
                        % Interpolate displacements
                        UNew(2*LengthOfCoordsOld+1) = 0.5*(U(2*(elementsFEM(refinementPatchtemp(2*tempi),3))-1)+U(2*(elementsFEM(refinementPatchtemp(2*tempi),1))-1));
                        UNew(2*LengthOfCoordsOld+2) = 0.5*(U(2*(elementsFEM(refinementPatchtemp(2*tempi),3)))+U(2*(elementsFEM(refinementPatchtemp(2*tempi),1))));                        
                        
                        FNew(4*LengthOfCoordsOld+1) = 0.5*(F(4*(elementsFEM(refinementPatchtemp(2*tempi),3))-3)+ ...
                                           F(4*(elementsFEM(refinementPatchtemp(2*tempi),1))-3));
                        FNew(4*LengthOfCoordsOld+2) = 0.5*(F(4*(elementsFEM(refinementPatchtemp(2*tempi),3))-2)+ ...
                                           F(4*(elementsFEM(refinementPatchtemp(2*tempi),1))-2));
                        FNew(4*LengthOfCoordsOld+3) = 0.5*(F(4*(elementsFEM(refinementPatchtemp(2*tempi),3))-1)+ ...
                                           F(4*(elementsFEM(refinementPatchtemp(2*tempi),1))-1));
                        FNew(4*LengthOfCoordsOld+4) = 0.5*(F(4*(elementsFEM(refinementPatchtemp(2*tempi),3)))+ ...
                                             F(4*(elementsFEM(refinementPatchtemp(2*tempi),1))));
                end
                
                
                switch longestEdgeNo
                    
                    case 1
                        elementsFEMNew(LengthOfElementsOld+1,:) = [elementsFEM(refinementPatchtemp(2*tempi-1),2:3), LengthOfCoordsOld+1];
                        elementsFEMNew(refinementPatch(tempi),:) = [elementsFEM(refinementPatchtemp(2*tempi-1),3), elementsFEM(refinementPatchtemp(2*tempi-1),1), LengthOfCoordsOld+1];
                       
                    case 2
                        elementsFEMNew(LengthOfElementsOld+1,:) = [elementsFEM(refinementPatchtemp(2*tempi-1),1:2), LengthOfCoordsOld+1];
                        elementsFEMNew(refinementPatch(tempi),:) = [elementsFEM(refinementPatchtemp(2*tempi-1),3), elementsFEM(refinementPatchtemp(2*tempi-1),1), LengthOfCoordsOld+1];
                        
                    case 3
                        elementsFEMNew(LengthOfElementsOld+1,:) = [elementsFEM(refinementPatchtemp(2*tempi-1),2:3), LengthOfCoordsOld+1];
                        elementsFEMNew(refinementPatch(tempi),:) = [elementsFEM(refinementPatchtemp(2*tempi-1),1:2), LengthOfCoordsOld+1];
                        
                end
                
                refinedEleIDList = [refinedEleIDList; refinementPatch(tempi); eleNeighbor4];
                
                
                
            else % Longest edge is on the boundary: elements can be directly refined itself, without any eleNeighbor4
                
                
                % Add the new added node onto dirichlet list
                dirichlet = LengthOfCoordsOld+1;
                
                switch longestEdgeNo
                    
                    case 1
                         
                        % Interpolate displacements
                        UNew(2*LengthOfCoordsOld+1) = 0.5*(U(2*(elementsFEM(refinementPatch(tempi),1))-1)+...
                                                            U(2*(elementsFEM(refinementPatch(tempi),2))-1));
                        UNew(2*LengthOfCoordsOld+2) = 0.5*(U(2*(elementsFEM(refinementPatch(tempi),1)))+ ...
                                                            U(2*(elementsFEM(refinementPatch(tempi),2))));
                         
                        FNew(4*LengthOfCoordsOld+1) = 0.5*(F(4*(elementsFEM(refinementPatch(tempi),1))-3)+...
                                                            F(4*(elementsFEM(refinementPatch(tempi),2))-3));
                        FNew(4*LengthOfCoordsOld+2) = 0.5*(F(4*(elementsFEM(refinementPatch(tempi),1))-2)+...
                                                            F(4*(elementsFEM(refinementPatch(tempi),2))-2));
                        FNew(4*LengthOfCoordsOld+3) = 0.5*(F(4*(elementsFEM(refinementPatch(tempi),1))-1)+...
                                                            F(4*(elementsFEM(refinementPatch(tempi),2))-1));
                        FNew(4*LengthOfCoordsOld+4) = 0.5*(F(4*(elementsFEM(refinementPatch(tempi),1)))+ ...
                                                            F(4*(elementsFEM(refinementPatch(tempi),2))));

                        elementsFEMNew(LengthOfElementsOld+1,:) = [elementsFEM(refinementPatch(tempi),2:3), LengthOfCoordsOld+1];
                        elementsFEMNew(refinementPatch(tempi),:) = [elementsFEM(refinementPatch(tempi),3), elementsFEM(refinementPatch(tempi),1), LengthOfCoordsOld+1];
                        
                        
                        
                    case 2
                         
                         % Interpolate displacements
                        UNew(2*LengthOfCoordsOld+1) = 0.5*(U(2*(elementsFEM(refinementPatch(tempi),2))-1)+U(2*(elementsFEM(refinementPatch(tempi),3))-1));
                        UNew(2*LengthOfCoordsOld+2) = 0.5*(U(2*(elementsFEM(refinementPatch(tempi),2)))+U(2*(elementsFEM(refinementPatch(tempi),3))));
                        
                        FNew(4*LengthOfCoordsOld+1) = 0.5*(F(4*(elementsFEM(refinementPatch(tempi),2))-3)+...
                                                            F(4*(elementsFEM(refinementPatch(tempi),3))-3));
                        FNew(4*LengthOfCoordsOld+2) = 0.5*(F(4*(elementsFEM(refinementPatch(tempi),2))-2)+...
                                                            F(4*(elementsFEM(refinementPatch(tempi),3))-2));
                        FNew(4*LengthOfCoordsOld+3) = 0.5*(F(4*(elementsFEM(refinementPatch(tempi),2))-1)+...
                                                            F(4*(elementsFEM(refinementPatch(tempi),3))-1));
                        FNew(4*LengthOfCoordsOld+4) = 0.5*(F(4*(elementsFEM(refinementPatch(tempi),2)))+ ...
                                                            F(4*(elementsFEM(refinementPatch(tempi),3))));

                        elementsFEMNew(LengthOfElementsOld+1,:) = [elementsFEM(refinementPatch(tempi),1:2), LengthOfCoordsOld+1];
                        elementsFEMNew(refinementPatch(tempi),:) = [elementsFEM(refinementPatch(tempi),3), elementsFEM(refinementPatch(tempi),1), LengthOfCoordsOld+1];
                        
                       
                    case 3
                         
                        % Interpolate displacements
                        UNew(2*LengthOfCoordsOld+1) = 0.5*(U(2*(elementsFEM(refinementPatch(tempi),3))-1)+U(2*(elementsFEM(refinementPatch(tempi),1))-1));
                        UNew(2*LengthOfCoordsOld+2) = 0.5*(U(2*(elementsFEM(refinementPatch(tempi),3)))+U(2*(elementsFEM(refinementPatch(tempi),1))));
                         
                        FNew(4*LengthOfCoordsOld+1) = 0.5*(F(4*(elementsFEM(refinementPatch(tempi),3))-3)+...
                                                        F(4*(elementsFEM(refinementPatch(tempi),1))-3));
                        FNew(4*LengthOfCoordsOld+2) = 0.5*(F(4*(elementsFEM(refinementPatch(tempi),3))-2)+...
                                                         F(4*(elementsFEM(refinementPatch(tempi),1))-2));
                        FNew(4*LengthOfCoordsOld+3) = 0.5*(F(4*(elementsFEM(refinementPatch(tempi),3))-1)+...
                                                         F(4*(elementsFEM(refinementPatch(tempi),1))-1));
                        FNew(4*LengthOfCoordsOld+4) = 0.5*(F(4*(elementsFEM(refinementPatch(tempi),3)))+ ...
                                                         F(4*(elementsFEM(refinementPatch(tempi),1))));

                        elementsFEMNew(LengthOfElementsOld+1,:) = [elementsFEM(refinementPatch(tempi),2:3), LengthOfCoordsOld+1];
                        elementsFEMNew(refinementPatch(tempi),:) = [elementsFEM(refinementPatch(tempi),1:2), LengthOfCoordsOld+1];
                        
                end
                
                refinedEleIDList = [refinedEleIDList;refinementPatch(tempi)];
                
                
            end
            
            
        end
        
    end
    
    
     
    
end

