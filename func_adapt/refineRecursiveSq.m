function [coordinatesFEMNew,elementsFEMNew,UNew,FNew,refinedEleIDList,dirichletNew,neumannNew,eleGenerationNew] = ...
    refineRecursiveSq(coordinatesFEM,elementsFEM,eleGeneration,U,F,refineT,winsize,FinishRefineRecursive,...
    CrackOrNot,CrackPath1,CrackPath2,CrackTip,CrackTipOrNot)


coordinatesFEMNew = coordinatesFEM; elementsFEMNew = elementsFEM; 
eleGenerationNew = eleGeneration; UNew = U; FNew = F; dirichletNew = 0; neumannNew = zeros(1,4);

refinedEleIDList = 0;
 
while FinishRefineRecursive ~= 0
    whileStopCriterion = 1; refineArray = refineT; ming = 1;
    
    while norm(whileStopCriterion) > 0
        
        SizeColOfRefineArraytemp2 = size(refineArray,2);
        refinementPatch = zeros(size(refineArray,1),1);
        
        for tempi = 1:size(refineArray,1)
            
            if whileStopCriterion(tempi) > 0
              
                refinementPatchtemp = getRefinementPatch(refineArray(tempi,SizeColOfRefineArraytemp2),eleGenerationNew,elementsFEMNew,coordinatesFEMNew,winsize);
                
                [sortg,indexg] = sort(eleGenerationNew(refinementPatchtemp));
                [rowtemp,~] = find(sortg==sortg(1));
                mingEle = refinementPatchtemp(indexg(rowtemp));
                 
                
                if  eleGenerationNew(refineArray(tempi,SizeColOfRefineArraytemp2)) == sortg(1)
                    
                    refinementPatch(tempi,1:size(refinementPatch,1)) = zeros(1,size(refinementPatch,1));
                    refinementPatch(tempi,1:length(mingEle)) = mingEle;
                     
                    whileStopCriterion(tempi) = 0;

                    ming(tempi) = sortg(1);

                      
                elseif eleGenerationNew(refineArray(tempi,SizeColOfRefineArraytemp2)) > sortg(1)
                    
                    SizeRowOfRefineArraytemp = size(refineArray,1);
                    SizeColOfRefineArraytemp = size(refineArray,2);
                    
                    tempRow = refineArray(tempi,1:SizeColOfRefineArraytemp);
                    refineArray(tempi,1:SizeColOfRefineArraytemp) = zeros(1,SizeColOfRefineArraytemp);
                     
                    % Delete all the zeros at the end
                    tempRowWithoutZeros = setdiff(tempRow,0,'stable');
                     
                    SizeColOftempRowWithoutZeros = size(tempRowWithoutZeros,2);
                    
                    refineArray(tempi,1:SizeColOftempRowWithoutZeros+1) = [tempRowWithoutZeros, mingEle(1)];
                    
                    whileStopCriterion(tempi) = 1;
                    
                    if length(mingEle) > 1
                        for tempj = 2:length(mingEle)
                             
                            refineArray(SizeRowOfRefineArraytemp+tempj-1,1:SizeColOftempRowWithoutZeros+1) = [tempRowWithoutZeros, mingEle(tempj)];
                            whileStopCriterion(SizeRowOfRefineArraytemp+tempj-1) = 1;
                        end
                        
                    end
                    
                end
                
            end
            
        end
        
        
    end
    
      
    % disp('======== INTO refinement part =========')
    
    coordinatesFEM = coordinatesFEMNew; elementsFEM = elementsFEMNew; eleGeneration = eleGenerationNew; U = UNew;
    
    % ====== Get refinement patch ========
    % ====== Compare ming and refine the min ming refinementPatch =======
     
    [sortming,indexming] = sort(ming);
     
    refinementPatchtemp = refinementPatch(indexming(1),:);
      
    
    for tempi = 1:length(refinementPatchtemp)
        
        if refinementPatchtemp(tempi) > 0
        
            coordx1 = coordinatesFEMNew(elementsFEMNew(refinementPatchtemp(tempi),1),1);
            coordy1 = coordinatesFEMNew(elementsFEMNew(refinementPatchtemp(tempi),1),2);
            coordx2 = coordinatesFEMNew(elementsFEMNew(refinementPatchtemp(tempi),2),1);
            coordy2 = coordinatesFEMNew(elementsFEMNew(refinementPatchtemp(tempi),2),2);
            coordx3 = coordinatesFEMNew(elementsFEMNew(refinementPatchtemp(tempi),3),1);
            coordy3 = coordinatesFEMNew(elementsFEMNew(refinementPatchtemp(tempi),3),2);
            coordx4 = coordinatesFEMNew(elementsFEMNew(refinementPatchtemp(tempi),4),1);
            coordy4 = coordinatesFEMNew(elementsFEMNew(refinementPatchtemp(tempi),4),2);
            
            % ------ Calculate nodes 5~8 coordinates -------
            coordx5 = 0.5*(coordx2+coordx3); coordy5 = 0.5*(coordy2+coordy3);
            coordx6 = 0.5*(coordx3+coordx4); coordy6 = 0.5*(coordy3+coordy4);
            coordx7 = 0.5*(coordx4+coordx1); coordy7 = 0.5*(coordy4+coordy1);
            coordx8 = 0.5*(coordx1+coordx2); coordy8 = 0.5*(coordy1+coordy2);
            coordx9 = 0.5*(coordx1+coordx3); coordy9 = 0.5*(coordy1+coordy3);
             
              
            % ------ node 5 -------
            % ------ Check points already exist or not ---------
            [row5,~] = find((coordinatesFEMNew(:,1) == coordx5) & (coordinatesFEMNew(:,2) == coordy5));
            
            if isempty(row5) == 0 % this mid point already exist; we don't need add new points,
            else % this mid point is a new point to generate
                LengthOfCoordsOld = size(coordinatesFEMNew,1);
                coordinatesFEMNew(LengthOfCoordsOld+1,1:2) = [coordx5,coordy5];
                row5 = LengthOfCoordsOld+1; 
                
                % !!! need extension for general cracks
                % UNew(2*row5-1) = nan; UNew(2*row5) = nan; FNew(4*row5-3:4*row5)=[0;0;0;0];
                % if CrackOrNot == 0 || coordx5-0.5*winsize > CrackTip(1) || coordx5+0.5*winsize < CrackTip(1) || coordy5-0.5*winsize > CrackTip(2)
                    UNew(2*row5-1) = 0.5 * ( UNew(2*elementsFEMNew(refinementPatchtemp(tempi),2)-1) + ...
                                            UNew(2*elementsFEMNew(refinementPatchtemp(tempi),3)-1) );
                    UNew(2*row5)   = 0.5 * ( UNew(2*elementsFEMNew(refinementPatchtemp(tempi),2)) + ...
                                            UNew(2*elementsFEMNew(refinementPatchtemp(tempi),3)) );
                    FNew(4*row5-3) = 0.5 * ( FNew(4*elementsFEMNew(refinementPatchtemp(tempi),2)-3) + ...
                                            FNew(4*elementsFEMNew(refinementPatchtemp(tempi),3)-3) );
                    FNew(4*row5-2) = 0.5 * ( FNew(4*elementsFEMNew(refinementPatchtemp(tempi),2)-2) + ...
                                            FNew(4*elementsFEMNew(refinementPatchtemp(tempi),3)-2) );
                    FNew(4*row5-1) = 0.5 * ( FNew(4*elementsFEMNew(refinementPatchtemp(tempi),2)-1) + ...
                                            FNew(4*elementsFEMNew(refinementPatchtemp(tempi),3)-1) );
                    FNew(4*row5-0) = 0.5 * ( FNew(4*elementsFEMNew(refinementPatchtemp(tempi),2)-0) + ...
                                            FNew(4*elementsFEMNew(refinementPatchtemp(tempi),3)-0) );
                % elseif coordx5 < CrackTip(1) % Left part is master
                     
                    
                
                % findEleNeighbors
                eleNeighbor = findRefineEleNeighbors(elementsFEMNew,refinementPatchtemp(tempi),5);
                  
                 
                if eleNeighbor == 0 % The new generated node is at the boundary
                    LengthOfDirichlettemp = length(dirichletNew);
                    dirichletNew(LengthOfDirichlettemp+1) = row5;
                    
                    LengthOfNeumanntemp = size(neumannNew,1);
                    neumannNew(LengthOfNeumanntemp+1,1:4) = [elementsFEMNew(refinementPatchtemp(tempi),2:3), row5, LengthOfDirichlettemp+1];
                    
                else % The new generated node is not at the boundary
                    % Modify the neighbor elements
                    elementsFEMNew(eleNeighbor,7) = row5;
                end
                
            end
              
            
            % ------ node 6 -------
            % ------ Check points already exist or not ---------
            [row6,~] = find((coordinatesFEMNew(:,1) == coordx6) & (coordinatesFEMNew(:,2) == coordy6));
            if isempty(row6) == 0 % this mid point already exist; we don't need add new points,
            else % this mid point is a new point to generate
                LengthOfCoordsOld = size(coordinatesFEMNew,1);
                coordinatesFEMNew(LengthOfCoordsOld+1,1:2) = [coordx6,coordy6];
                row6 = LengthOfCoordsOld+1;
                
                % UNew(2*row6-1) = nan; UNew(2*row6) = nan; FNew(4*row6-3:4*row6)=[0;0;0;0];
                
                % !!! need extension for general cracks
               
                if CrackOrNot == 0 || coordx6-0.5*winsize > CrackTip(1) || coordx6+0.5*winsize < CrackTip(1) || coordy6-0.5*winsize > CrackTip(2)
                
                    UNew(2*row6-1) = 0.5 * ( UNew(2*elementsFEMNew(refinementPatchtemp(tempi),3)-1) + ...
                                             UNew(2*elementsFEMNew(refinementPatchtemp(tempi),4)-1) );
                    UNew(2*row6)   = 0.5 * ( UNew(2*elementsFEMNew(refinementPatchtemp(tempi),3)) + ...
                                             UNew(2*elementsFEMNew(refinementPatchtemp(tempi),4)) );

                    FNew(4*row6-3) = 0.5 * ( FNew(4*elementsFEMNew(refinementPatchtemp(tempi),3)-3) + ...
                                            FNew(4*elementsFEMNew(refinementPatchtemp(tempi),4)-3) );
                    FNew(4*row6-2) = 0.5 * ( FNew(4*elementsFEMNew(refinementPatchtemp(tempi),3)-2) + ...
                                            FNew(4*elementsFEMNew(refinementPatchtemp(tempi),4)-2) );
                    FNew(4*row6-1) = 0.5 * ( FNew(4*elementsFEMNew(refinementPatchtemp(tempi),3)-1) + ...
                                            FNew(4*elementsFEMNew(refinementPatchtemp(tempi),4)-1) );
                    FNew(4*row6-0) = 0.5 * ( FNew(4*elementsFEMNew(refinementPatchtemp(tempi),3)-0) + ...
                                            FNew(4*elementsFEMNew(refinementPatchtemp(tempi),4)-0) );
                                        
                elseif coordx6 < CrackTip(1) % Left part is master, and right part will be discarded
                    UNew(2*row6-1) = UNew(2*elementsFEMNew(refinementPatchtemp(tempi),4)-1) ;
                    UNew(2*row6)   = UNew(2*elementsFEMNew(refinementPatchtemp(tempi),4)) ;
                    FNew(4*row6-3) = FNew(4*elementsFEMNew(refinementPatchtemp(tempi),4)-3) ;
                    FNew(4*row6-2) = FNew(4*elementsFEMNew(refinementPatchtemp(tempi),4)-2) ;
                    FNew(4*row6-1) = FNew(4*elementsFEMNew(refinementPatchtemp(tempi),4)-1) ;
                    FNew(4*row6-0) = FNew(4*elementsFEMNew(refinementPatchtemp(tempi),4)-0) ;
                elseif coordx6 > CrackTip(1) % Right part is master, and left part will be discarded
                    UNew(2*row6-1) = UNew(2*elementsFEMNew(refinementPatchtemp(tempi),3)-1) ;
                    UNew(2*row6)   = UNew(2*elementsFEMNew(refinementPatchtemp(tempi),3)) ;
                    FNew(4*row6-3) = FNew(4*elementsFEMNew(refinementPatchtemp(tempi),3)-3) ;
                    FNew(4*row6-2) = FNew(4*elementsFEMNew(refinementPatchtemp(tempi),3)-2) ;
                    FNew(4*row6-1) = FNew(4*elementsFEMNew(refinementPatchtemp(tempi),3)-1) ;
                    FNew(4*row6-0) = FNew(4*elementsFEMNew(refinementPatchtemp(tempi),3)-0) ;
                end
                    
                    
                    
                 
                % findEleNeighbors
                eleNeighbor = findRefineEleNeighbors(elementsFEMNew,refinementPatchtemp(tempi),6);
                 
                if eleNeighbor == 0 % The new generated node is at the boundary
                    LengthOfDirichlettemp = length(dirichletNew);
                    dirichletNew(LengthOfDirichlettemp+1) = row6;
                    
                    LengthOfNeumanntemp = size(neumannNew,1);
                    neumannNew(LengthOfNeumanntemp+1,1:4) = [elementsFEMNew(refinementPatchtemp(tempi),3:4), row6, LengthOfDirichlettemp+1];
                     
                else % The new generated node is not at the boundary
                    % Modify the neighbor elements
                    elementsFEMNew(eleNeighbor,8) = row6;
                     
                end
            end
            
            
            % ------ node 7 -------
            % ------ Check points already exist or not ---------
            [row7,~] = find((coordinatesFEMNew(:,1) == coordx7) & (coordinatesFEMNew(:,2) == coordy7));
            if isempty(row7) == 0 % this mid point already exist; we don't need add new points,
            else % this mid point is a new point to generate
                LengthOfCoordsOld = size(coordinatesFEMNew,1);
                coordinatesFEMNew(LengthOfCoordsOld+1,1:2) = [coordx7,coordy7];
                row7 = LengthOfCoordsOld+1;
                
                % UNew(2*row7-1) = nan; UNew(2*row7) = nan; FNew(4*row7-3:4*row7)=[0;0;0;0];
                
                UNew(2*row7-1) = 0.5 * ( UNew(2*elementsFEMNew(refinementPatchtemp(tempi),4)-1) + ...
                                         UNew(2*elementsFEMNew(refinementPatchtemp(tempi),1)-1) );
                UNew(2*row7)   = 0.5 * ( UNew(2*elementsFEMNew(refinementPatchtemp(tempi),4)) + ...
                                         UNew(2*elementsFEMNew(refinementPatchtemp(tempi),1)) );
                FNew(4*row7-3) = 0.5 * ( FNew(4*elementsFEMNew(refinementPatchtemp(tempi),4)-3) + ...
                                        FNew(4*elementsFEMNew(refinementPatchtemp(tempi),1)-3) );
                FNew(4*row7-2) = 0.5 * ( FNew(4*elementsFEMNew(refinementPatchtemp(tempi),4)-2) + ...
                                        FNew(4*elementsFEMNew(refinementPatchtemp(tempi),1)-2) );
                FNew(4*row7-1) = 0.5 * ( FNew(4*elementsFEMNew(refinementPatchtemp(tempi),4)-1) + ...
                                        FNew(4*elementsFEMNew(refinementPatchtemp(tempi),1)-1) );
                FNew(4*row7-0) = 0.5 * ( FNew(4*elementsFEMNew(refinementPatchtemp(tempi),4)-0) + ...
                                        FNew(4*elementsFEMNew(refinementPatchtemp(tempi),1)-0) );
                                    
                % findEleNeighbors
                eleNeighbor = findRefineEleNeighbors(elementsFEMNew,refinementPatchtemp(tempi),7);
                   
                  
                if eleNeighbor == 0 % The new generated node is at the boundary
                    LengthOfDirichlettemp = length(dirichletNew);
                    dirichletNew(LengthOfDirichlettemp+1) = row7;
                    
                    LengthOfNeumanntemp = size(neumannNew,1);
                    neumannNew(LengthOfNeumanntemp+1,1:4) = [elementsFEMNew(refinementPatchtemp(tempi),[4,1]), row7, LengthOfDirichlettemp+1];
                    
                else % The new generated node is not at the boundary
                    % Modify the neighbor elements
                    elementsFEMNew(eleNeighbor,5) = row7;
                      
                end
            end
            
            
            % ------ node 8 -------
            % ------ Check points already exist or not ---------
            [row8,~] = find((coordinatesFEMNew(:,1) == coordx8) & (coordinatesFEMNew(:,2) == coordy8));
            if isempty(row8) == 0 % this mid point already exist; we don't need add new points,
            else % this mid point is a new point to generate
                LengthOfCoordsOld = size(coordinatesFEMNew,1);
                coordinatesFEMNew(LengthOfCoordsOld+1,1:2) = [coordx8,coordy8];
                row8 = LengthOfCoordsOld+1;
                
                % UNew(2*row8-1) = nan; UNew(2*row8) = nan; FNew(4*row8-3:4*row8)=[0;0;0;0];
                
                % !!! need extension for general cracks
                if CrackOrNot == 0 || coordx8-0.5*winsize > CrackTip(1) || coordx8+0.5*winsize < CrackTip(1) || coordy8-0.5*winsize > CrackTip(2)
                
                    UNew(2*row8-1) = 0.5 * ( UNew(2*elementsFEMNew(refinementPatchtemp(tempi),1)-1) + ...
                        UNew(2*elementsFEMNew(refinementPatchtemp(tempi),2)-1) );
                    UNew(2*row8)   = 0.5 * ( UNew(2*elementsFEMNew(refinementPatchtemp(tempi),1)) + ...
                        UNew(2*elementsFEMNew(refinementPatchtemp(tempi),2)) );
                    FNew(4*row8-3) = 0.5 * ( FNew(4*elementsFEMNew(refinementPatchtemp(tempi),1)-3) + ...
                        FNew(4*elementsFEMNew(refinementPatchtemp(tempi),2)-3) );
                    FNew(4*row8-2) = 0.5 * ( FNew(4*elementsFEMNew(refinementPatchtemp(tempi),1)-2) + ...
                        FNew(4*elementsFEMNew(refinementPatchtemp(tempi),2)-2) );
                    FNew(4*row8-1) = 0.5 * ( FNew(4*elementsFEMNew(refinementPatchtemp(tempi),1)-1) + ...
                        FNew(4*elementsFEMNew(refinementPatchtemp(tempi),2)-1) );
                    FNew(4*row8-0) = 0.5 * ( FNew(4*elementsFEMNew(refinementPatchtemp(tempi),1)-0) + ...
                        FNew(4*elementsFEMNew(refinementPatchtemp(tempi),2)-0) );
                                        
                elseif coordx8 < CrackTip(1) % Left part is master, and right part will be discarded
                    UNew(2*row8-1) = UNew(2*elementsFEMNew(refinementPatchtemp(tempi),1)-1) ;
                    UNew(2*row8)   = UNew(2*elementsFEMNew(refinementPatchtemp(tempi),1)) ;
                    FNew(4*row8-3) = FNew(4*elementsFEMNew(refinementPatchtemp(tempi),1)-3) ;
                    FNew(4*row8-2) = FNew(4*elementsFEMNew(refinementPatchtemp(tempi),1)-2) ;
                    FNew(4*row8-1) = FNew(4*elementsFEMNew(refinementPatchtemp(tempi),1)-1) ;
                    FNew(4*row8-0) = FNew(4*elementsFEMNew(refinementPatchtemp(tempi),1)-0) ;
                elseif coordx8 > CrackTip(1) % Right part is master, and left part will be discarded
                    UNew(2*row8-1) = UNew(2*elementsFEMNew(refinementPatchtemp(tempi),2)-1) ;
                    UNew(2*row8)   = UNew(2*elementsFEMNew(refinementPatchtemp(tempi),2)) ;
                    FNew(4*row8-3) = FNew(4*elementsFEMNew(refinementPatchtemp(tempi),2)-3) ;
                    FNew(4*row8-2) = FNew(4*elementsFEMNew(refinementPatchtemp(tempi),2)-2) ;
                    FNew(4*row8-1) = FNew(4*elementsFEMNew(refinementPatchtemp(tempi),2)-1) ;
                    FNew(4*row8-0) = FNew(4*elementsFEMNew(refinementPatchtemp(tempi),2)-0) ;
                end
                
                
                                     
                % findEleNeighbors
                eleNeighbor = findRefineEleNeighbors(elementsFEMNew,refinementPatchtemp(tempi),8);
                 
                
                if eleNeighbor == 0 % The new generated node is at the boundary
                    LengthOfDirichlettemp = length(dirichletNew);
                    dirichletNew(LengthOfDirichlettemp+1) = row8;
                    
                    LengthOfNeumanntemp = size(neumannNew,1);
                    neumannNew(LengthOfNeumanntemp+1,1:4) = [elementsFEMNew(refinementPatchtemp(tempi),1:2), row8, LengthOfDirichlettemp+1];
                    
                else % The new generated node is not at the boundary
                    % Modify the neighbor elements
                    elementsFEMNew(eleNeighbor,6) = row8;
                     
                end
            end
              
            
            % ------ node 9 -------
            LengthOfCoordsOld = size(coordinatesFEMNew,1);
            coordinatesFEMNew(LengthOfCoordsOld+1,1:2) = [coordx9,coordy9];
            row9 = LengthOfCoordsOld+1;
            
            % UNew(2*row9-1) = nan; UNew(2*row9) = nan; FNew(4*row9-3:4*row9)=[0;0;0;0];
            % !!! need extension for general cracks
                if CrackOrNot == 0 || coordx9-0.5*winsize > CrackTip(1) || coordx9+0.5*winsize < CrackTip(1) || coordy9-0.5*winsize > CrackTip(2)
                
                    UNew(2*row9-1) = 0.25 * ( UNew(2*elementsFEMNew(refinementPatchtemp(tempi),1)-1) + ...
                        UNew(2*elementsFEMNew(refinementPatchtemp(tempi),2)-1) + ...
                        UNew(2*elementsFEMNew(refinementPatchtemp(tempi),3)-1) + ...
                        UNew(2*elementsFEMNew(refinementPatchtemp(tempi),4)-1) );
                    UNew(2*row9)   = 0.25 * ( UNew(2*elementsFEMNew(refinementPatchtemp(tempi),1)) + ...
                        UNew(2*elementsFEMNew(refinementPatchtemp(tempi),2)) + ...
                        UNew(2*elementsFEMNew(refinementPatchtemp(tempi),3)) + ...
                        UNew(2*elementsFEMNew(refinementPatchtemp(tempi),4)) );
                    FNew(4*row9-3) = 0.25 * ( FNew(4*elementsFEMNew(refinementPatchtemp(tempi),1)-3) + ...
                        FNew(4*elementsFEMNew(refinementPatchtemp(tempi),2)-3) + ...
                        FNew(4*elementsFEMNew(refinementPatchtemp(tempi),3)-3)+ ...
                        FNew(4*elementsFEMNew(refinementPatchtemp(tempi),4)-3));
                    FNew(4*row9-2) = 0.25 * ( FNew(4*elementsFEMNew(refinementPatchtemp(tempi),1)-2) + ...
                        FNew(4*elementsFEMNew(refinementPatchtemp(tempi),2)-2) + ...
                        FNew(4*elementsFEMNew(refinementPatchtemp(tempi),3)-2)+ ...
                        FNew(4*elementsFEMNew(refinementPatchtemp(tempi),4)-2));
                    FNew(4*row9-1) = 0.25 * ( FNew(4*elementsFEMNew(refinementPatchtemp(tempi),1)-1) + ...
                        FNew(4*elementsFEMNew(refinementPatchtemp(tempi),2)-1) + ...
                        FNew(4*elementsFEMNew(refinementPatchtemp(tempi),3)-1)+ ...
                        FNew(4*elementsFEMNew(refinementPatchtemp(tempi),4)-1));
                    FNew(4*row9-0) = 0.25 * ( FNew(4*elementsFEMNew(refinementPatchtemp(tempi),1)-0) + ...
                        FNew(4*elementsFEMNew(refinementPatchtemp(tempi),2)-0) + ...
                        FNew(4*elementsFEMNew(refinementPatchtemp(tempi),3)-0)+ ...
                        FNew(4*elementsFEMNew(refinementPatchtemp(tempi),4)-0));
                                        
                elseif coordx9 < CrackTip(1) % Left part is master, and right part will be discarded
                    UNew(2*row9-1) = 0.5 * ( UNew(2*elementsFEMNew(refinementPatchtemp(tempi),1)-1) + ...
                        UNew(2*elementsFEMNew(refinementPatchtemp(tempi),4)-1) );
                    UNew(2*row9)   = 0.5 * ( UNew(2*elementsFEMNew(refinementPatchtemp(tempi),1)) + ...
                        UNew(2*elementsFEMNew(refinementPatchtemp(tempi),4)) );
                    FNew(4*row9-3) = 0.5 * ( FNew(4*elementsFEMNew(refinementPatchtemp(tempi),1)-3) + ...
                        FNew(4*elementsFEMNew(refinementPatchtemp(tempi),4)-3));
                    FNew(4*row9-2) = 0.5 * ( FNew(4*elementsFEMNew(refinementPatchtemp(tempi),1)-2) + ...
                        FNew(4*elementsFEMNew(refinementPatchtemp(tempi),4)-2));
                    FNew(4*row9-1) = 0.5 * ( FNew(4*elementsFEMNew(refinementPatchtemp(tempi),1)-1) + ...
                        FNew(4*elementsFEMNew(refinementPatchtemp(tempi),4)-1));
                    FNew(4*row9-0) = 0.5 * ( FNew(4*elementsFEMNew(refinementPatchtemp(tempi),1)-0) + ...
                        FNew(4*elementsFEMNew(refinementPatchtemp(tempi),4)-0));
                    
                elseif coordx9 > CrackTip(1) % Right part is master, and left part will be discarded
                    UNew(2*row9-1) = 0.5 * ( UNew(2*elementsFEMNew(refinementPatchtemp(tempi),2)-1) + ...
                        UNew(2*elementsFEMNew(refinementPatchtemp(tempi),3)-1) );
                    UNew(2*row9)   = 0.5 * ( UNew(2*elementsFEMNew(refinementPatchtemp(tempi),2)) + ...
                        UNew(2*elementsFEMNew(refinementPatchtemp(tempi),3)) );
                    FNew(4*row9-3) = 0.5 * ( FNew(4*elementsFEMNew(refinementPatchtemp(tempi),2)-3) + ...
                        FNew(4*elementsFEMNew(refinementPatchtemp(tempi),3)-3));
                    FNew(4*row9-2) = 0.5 * ( FNew(4*elementsFEMNew(refinementPatchtemp(tempi),2)-2) + ...
                        FNew(4*elementsFEMNew(refinementPatchtemp(tempi),3)-2));
                    FNew(4*row9-1) = 0.5 * ( FNew(4*elementsFEMNew(refinementPatchtemp(tempi),2)-1) + ...
                        FNew(4*elementsFEMNew(refinementPatchtemp(tempi),3)-1));
                    FNew(4*row9-0) = 0.5 * ( FNew(4*elementsFEMNew(refinementPatchtemp(tempi),2)-0) + ...
                        FNew(4*elementsFEMNew(refinementPatchtemp(tempi),3)-0));
                    
                end
            
            
            
            
            % Refine elements
             
            LengthOfElementsOld = size(elementsFEMNew,1);
            elementsFEMNew(LengthOfElementsOld+1,:) = [row8,elementsFEM(refinementPatchtemp(tempi),2),row5,row9,0,0,0,0];
            elementsFEMNew(LengthOfElementsOld+2,:) = [row9,row5,elementsFEM(refinementPatchtemp(tempi),3),row6,0,0,0,0];
            elementsFEMNew(LengthOfElementsOld+3,:) = [row7,row9,row6,elementsFEM(refinementPatchtemp(tempi),4),0,0,0,0];
            elementsFEMNew(refinementPatchtemp(tempi),:) = [elementsFEM(refinementPatchtemp(tempi),1), row8, row9, row7,0,0,0,0];
            
            eleGenerationNew(LengthOfElementsOld+1) = eleGeneration(refinementPatchtemp(tempi))+1;
            eleGenerationNew(LengthOfElementsOld+2) = eleGeneration(refinementPatchtemp(tempi))+1;
            eleGenerationNew(LengthOfElementsOld+3) = eleGeneration(refinementPatchtemp(tempi))+1;
            eleGenerationNew(refinementPatchtemp(tempi),:) = eleGeneration(refinementPatchtemp(tempi))+1;
            
            lengthOfRefinedEleIDList = length(refinedEleIDList);
            refinedEleIDList(lengthOfRefinedEleIDList+1) =  refinementPatchtemp(tempi) ;
     
        end
        
    end
     
   
    
    FinishRefineRecursive = 1-ismember(refineT,refinementPatchtemp);
    
    
end

refinedEleIDList = refinedEleIDList';











