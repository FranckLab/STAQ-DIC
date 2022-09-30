%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function AL-DIC pre-Subproblem 2
% Object: to find N_{cut} and N_{tip}
% Author: Jin Yang
% Last date modified: 2018.03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [EnrHAndTipEleIndex,EnrTipEleIndex] = EnrEleWithFrac(coordinates,elements,CrackPath1,CrackPath2,CrackTip)

% ------- Info of CrackTip path-line --------
k_tip = -0.5*(CrackPath1(1)/CrackPath1(2)+CrackPath2(1)/CrackPath2(2));
CrackPathCen = [k_tip/(CrackTip(2)-k_tip*CrackTip(1)), -1/(CrackTip(2)-k_tip*CrackTip(1))];

% -------------------------------------------
EnrHAndTipEleIndex = []; EnrTipEleIndex = [];
% -------------------------------------------
elements = [elements,zeros(size(elements,1),1)];
% Add additional column onto elements, to indicate normal elements or
% enriched elements with jump-only or jump+CrackMode sigular functions.
% [elements(1,:),   0;
%  elements(2,:),   0;
%  elements(3,:),   0;
%  ......
%  elements(end,:), 0];
% After modification, the last column means
% 0 : Normal element without crack passing through
% 1 : With crack passing through but not contain crack tip
% 2 : Include crack tip 
% -------------------------------------------

for tempi = 1:size(elements,1)
    
    % ------------------------------
    % Find each element region
    x = coordinates(elements(tempi,1:4),1); x0 = 0.5*(x(1)+x(3));
    y = coordinates(elements(tempi,1:4),2); y0 = 0.5*(y(1)+y(3));
    winsize = x(3)-x(1); % If not square, winsizey = y(3)-y(1);
    
    % ------------------------------
    % Crack pass each element or not
    PassCrackOrNot = 0; % Initialize as NOT passing the local subset
    
    % if y0-winsize*0.5 < CrackTip(2) % Need to modify for general crack paths
        tempCen_1 = CrackPathCen(1)*(x0-winsize/2) + CrackPathCen(2)*(y0-winsize/2) + 1; % Node 1
        tempCen_3 = CrackPathCen(1)*(x0+winsize/2) + CrackPathCen(2)*(y0+winsize/2) + 1; % Node 3
        tempCen_2 = CrackPathCen(1)*(x0+winsize/2) + CrackPathCen(2)*(y0-winsize/2) + 1; % Node 2
        tempCen_4 = CrackPathCen(1)*(x0-winsize/2) + CrackPathCen(2)*(y0+winsize/2) + 1; % Node 4
        
        temp1_1 = CrackPath1(1)*(x0-winsize/2) + CrackPath1(2)*(y0-winsize/2) + 1; % Node 1
        temp1_3 = CrackPath1(1)*(x0+winsize/2) + CrackPath1(2)*(y0+winsize/2) + 1; % Node 3
        temp1_2 = CrackPath1(1)*(x0+winsize/2) + CrackPath1(2)*(y0-winsize/2) + 1; % Node 2
        temp1_4 = CrackPath1(1)*(x0-winsize/2) + CrackPath1(2)*(y0+winsize/2) + 1; % Node 4
        
        temp2_1 = CrackPath2(1)*(x0-winsize/2) + CrackPath2(2)*(y0-winsize/2) + 1; % Node 1
        temp2_3 = CrackPath2(1)*(x0+winsize/2) + CrackPath2(2)*(y0+winsize/2) + 1; % Node 3
        temp2_2 = CrackPath2(1)*(x0+winsize/2) + CrackPath2(2)*(y0-winsize/2) + 1; % Node 2
        temp2_4 = CrackPath2(1)*(x0-winsize/2) + CrackPath2(2)*(y0+winsize/2) + 1; % Node 4

        if ((tempCen_1*tempCen_3 < 0)||(tempCen_2*tempCen_4) < 0) ...
                || ((temp1_1*temp1_3 < 0)||(temp1_2*temp1_4) < 0) || ((temp2_1*temp2_3 < 0)||(temp2_2*temp2_4) < 0)
            elements(tempi,end) = 1; % PassCrackOrNot = 1;
            if y0-winsize*0.5 < CrackTip(2) && y0+winsize*0.5 > CrackTip(2) && x0-winsize*0.5 < CrackTip(1) && x0+winsize*0.5 > CrackTip(1)
                elements(tempi,end) = 2; % PassCrackOrNot = 1;
            end
            
        end
    % end
    
end

% ------ Collect all the enriched elements ------
[row,~] = find(elements(:,end)==1); EnrHNodes = unique(elements(row,1:8)); 
[row,~] = find(elements(:,end)==2); EnrTipNodes = unique(elements(row,1:8));
EnrHAndTipNodes = union(EnrHNodes,EnrTipNodes);

% ------ Assign enriched H elements nodes No ------
lengthOfCoords = size(coordinates,1); EnrEleIndexNo = 1;
for tempk = 1:length(EnrHAndTipNodes)
    if EnrHAndTipNodes(tempk) ~= 0
        EnrHAndTipEleIndex(EnrEleIndexNo, 1:2) = [EnrHAndTipNodes(tempk), lengthOfCoords+1];
        lengthOfCoords = lengthOfCoords + 1; EnrEleIndexNo = EnrEleIndexNo + 1;
    end
end
 
% ------ Assign enriched tip elements nodes No ------
EnrEleIndexNo = 1;
for tempk = 1:length(EnrTipNodes)
    if EnrTipNodes(tempk) ~= 0
        EnrTipEleIndex(EnrEleIndexNo, 1:2) = [EnrTipNodes(tempk), lengthOfCoords+1];
        lengthOfCoords = lengthOfCoords + 4;  % Plus 4 is for four NC functions
        EnrEleIndexNo = EnrEleIndexNo + 1;
    end
end



 