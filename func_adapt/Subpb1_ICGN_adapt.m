%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function AL-DIC Subproblem 1
% Object: to find deformation field using local methods
% Author: Jin Yang
% Last date modified: 2018.03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [U, H, stepwithinwhile, PassCrackOrNot] = ...
        Subpb1_ICGN( x0, y0, f, g, DfDx, DfDy, DfAxis, winsize,...
                H, beta, mu, udual, vdual, UOld, FOld, tol, method, ...
                CrackOrNot,CrackPath1,CrackPath2,CrackTip)
 
DfDxStartx = DfAxis(1); DfDxStarty = DfAxis(3);
            
% ====== Info of CrackTip path-line ======
k_tip = -0.5*(CrackPath1(1)/CrackPath1(2)+CrackPath2(1)/CrackPath2(2));
CrackPathCen = [k_tip/(CrackTip(2)-k_tip*CrackTip(1)), -1/(CrackTip(2)-k_tip*CrackTip(1))];

% ---------------------------
% Find local subset region
x  = [x0-winsize/2 ; x0+winsize/2 ; x0+winsize/2 ; x0-winsize/2]; % [coordinates(elements(j,:),1)];
y  = [y0-winsize/2 ; y0-winsize/2 ; y0+winsize/2 ; y0+winsize/2];  % [coordinates(elements(j,:),2)];

% ---------------------------
% Crack pass this subset or not
PassCrackOrNot = 0; % Initialize as NOT passing the local subset
if CrackOrNot == 1
    tempCen = CrackPathCen*[x';y']+ ones(1,4);  
    if (tempCen(1)*tempCen(3) < 0) || (tempCen(2)*tempCen(4) < 0) 
        PassCrackOrNot = 1; % If Pass cracks, finding the master part; Unfinished here
    end
end
    
% ---------------------------
% Initialization: Get P0
P0 = [FOld(1) FOld(2) FOld(3) FOld(4) UOld(1) UOld(2)]'; P = P0;  
 
% ---------------------------
% Find region for f
tempf = zeros(winsize+1,winsize+1);
tempf(1:(winsize+1),1:(winsize+1)) = f(x(1):x(3), y(1):y(3));
if norm(H) < abs(eps) % In the first step, H is zeros
    H = zeros(6,6);
    [tempCoordy, tempCoordx] = meshgrid(y(1):y(3),x(1):x(3));
    tempCoordx = tempCoordx(:); tempCoordy = tempCoordy(:);
    
    % Compute for all the pixels in the master part
    for tempij = 1:size(tempCoordx,1)
        if PassCrackOrNot == 1 % Crack pass through this local subset
            % if x0 < CrackTip(1) % Left part is master, and right part will be discarded
            if CrackPathCen(1)*x0 + CrackPathCen(2)*y0 + 1 > 0 % y0 < CrackTip(2) % Bottom part is master, and top part will be discarded
                % if tempCoordx(tempij) > CrackTip(1) && tempCoordy(tempij) < CrackTip(2)
                %     tempCoordx(tempij) = tempCoordx(tempij)-winsize;
                % end
                % if tempCoordx(tempij)+1-x(1)>1 && tempCoordx(tempij)>x0-0.8*winsize
                %     H = H + ([DfDx(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty) DfDy(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty)]*...
                %         [tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1 0; 0 tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1])'* ...
                %         ([DfDx(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty) DfDy(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty)]*...
                %         [tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1 0; 0 tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1]);
                % end
                if tempCoordx(tempij)*CrackPathCen(1) + tempCoordy(tempij)*CrackPathCen(2) + 1 < 0
                    tempCoordy(tempij) = tempCoordy(tempij)-winsize;
                end
                if (tempCoordy(tempij)+1-y(1) > 1) && (tempCoordy(tempij)+1-y(1) < size(tempf,2)) ...
                        && (tempCoordx(tempij)*CrackPath2(1) + tempCoordy(tempij)*CrackPath2(2) + 1 > 0)
                    H = H + ([DfDx(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty) DfDy(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty)]*...
                        [tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1 0; 0 tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1])'* ...
                        ([DfDx(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty) DfDy(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty)]*...
                        [tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1 0; 0 tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1]);
                end
            % elseif (x0 >= CrackTip(1)) % Right part is master, and left part will be discarded
            elseif CrackPathCen(1)*x0 + CrackPathCen(2)*y0 + 1 <= 0 % Top part is master, and bottom part will be discarded
                % if tempCoordx(tempij) < CrackTip(1) && tempCoordy(tempij) < CrackTip(2)
                %     tempCoordx(tempij) = tempCoordx(tempij)+winsize;
                % end
                % if tempCoordx(tempij)<size(DfDx,1) && tempCoordx(tempij)+1-x(1)<size(tempf,1) && tempCoordx(tempij)<x0+0.5*winsize
                %     H = H + ([DfDx(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty) DfDy(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty)]*...
                %         [tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1 0; 0 tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1])'* ...
                %         ([DfDx(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty) DfDy(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty)]*...
                %         [tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1 0; 0 tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1]);
                % end
                if tempCoordx(tempij)*CrackPathCen(1) + tempCoordy(tempij)*CrackPathCen(2) + 1 > 0
                    tempCoordy(tempij) = tempCoordy(tempij)+winsize;
                end
                if (tempCoordy(tempij)-DfDxStarty < size(DfDy,2)) && (tempCoordy(tempij)+1-y(1) < size(tempf,2)) ...
                        && (tempCoordx(tempij)*CrackPath1(1) + tempCoordy(tempij)*CrackPath1(2) + 1 < 0)
                    H = H + ([DfDx(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty) DfDy(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty)]*...
                        [tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1 0; 0 tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1])'* ...
                        ([DfDx(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty) DfDy(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty)]*...
                        [tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1 0; 0 tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1]);
                end
            end
        else % Crack didn't pass through this local subset
            H = H + ([DfDx(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty) DfDy(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty)]*...
                [tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1 0; 0 tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1])'* ...
                ([DfDx(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty) DfDy(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty)]*...
                [tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1 0; 0 tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1]);
        end
    end
end

meanf = mean(tempf(:)); bottomf = sqrt((length(tempf(:))-1)*var(tempf(:)));
H2 = H(5:6,5:6)*2/(bottomf^2) + [mu 0; 0 mu];
%%%%%%%%%%%%%%%% Tried below, not succeed %%%%%%%%%%%%%%%%%%%
% if CrackOrNot == 0
%     H2 = H(5:6,5:6)*2/(bottomf^2) + [mu 0; 0 mu];
% else
%     H2 = H*2/(bottomf^2) + [beta 0 0 0 0 0; 0 beta 0 0 0 0; 0 0 beta 0 0 0; 0 0 0 beta 0 0 ; 0 0 0 0 mu 0; 0 0 0 0 0 mu];
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --------------------------
% Initialize while loop
normOfWOld = 2; normOfWNew = 1; stepwithinwhile=0;
switch method   % For Gauss-Newton method
    case 'GaussNewton'
        delta = 0;
    case 'LevenbergMarquardt'
        delta = 0.001; % For Levenberg-Marquardt method 
        KappaOld = 1e6; KappaNew = 1e6; KappaStore = zeros(10,1); PStore = zeros(10,6);
    otherwise
        delta = 0;
end

while(  stepwithinwhile <= 30 && normOfWNew > tol  )
           
    stepwithinwhile = stepwithinwhile+1;

    % Find region for g
    [tempCoordy, tempCoordx] = meshgrid(y(1):y(3),x(1):x(3));
 
    tempCoordx = tempCoordx - x0*ones(winsize+1,winsize+1);
    tempCoordy = tempCoordy - y0*ones(winsize+1,winsize+1);
    u22 = (1+P(1))*tempCoordx + P(3)*tempCoordy + (x0+P(5))*ones(winsize+1,winsize+1);
    v22 = P(2)*tempCoordx + (1+P(4))*tempCoordy + (y0+P(6))*ones(winsize+1,winsize+1);
   
    row1 = find(u22<3); row2 = find(u22>size(f,1)-2); row3 = find(v22<3); row4 = find(v22>size(f,2)-2); 
    if ~isempty([row1; row2; row3; row4])
        normOfWNew = 1e6;
        warning('Out of image boundary!')
        break;
    else
    
        tempg = zeros(size(tempf,1)*size(tempf,2),1);
        [tempCoordy, tempCoordx] = meshgrid(1:winsize+1,1:winsize+1);
        tempCoordx = tempCoordx(:); tempCoordy = tempCoordy(:);
  
        for tempij = 1:size(tempCoordx,1)
            tempg(tempij,1)= ...
                fungInterpolation_g(  u22(tempCoordx(tempij),tempCoordy(tempij)) , v22(tempCoordx(tempij),tempCoordy(tempij)), ...
                g((floor(u22(tempCoordx(tempij),tempCoordy(tempij)))-1):(floor(u22(tempCoordx(tempij),tempCoordy(tempij)))+2), ...
                (floor(v22(tempCoordx(tempij),tempCoordy(tempij)))-1):(floor(v22(tempCoordx(tempij),tempCoordy(tempij)))+2))  );
        end
        
        tempg = reshape(tempg, winsize+1, winsize+1);
        meang = mean(tempg(:));
        bottomg = sqrt((length(tempg(:))-1)*var(tempg(:)));

        % ============ For Levenberg-Marquardt method ============    
        switch method
            case 'LevenbergMarquardt'
            % Compute functinoal error
            KappaOld = KappaNew;
            Kappatemp = (tempf-meanf)/bottomf - (tempg-meang)/bottomg;
            Kappatemp = Kappatemp.*Kappatemp;
            KappaNew = sum(Kappatemp(:));  

               if KappaNew < 1.02*KappaOld
                   delta = delta/10; 
               else
                   delta = delta*10;
                   % Perform P inverse 
                   DeltaP = -DeltaP;
                   tempP1 =  (-DeltaP(1)-DeltaP(1)*DeltaP(4)+DeltaP(2)*DeltaP(3))/tempDelP;
                   tempP2 =  -DeltaP(2)/tempDelP;
                   tempP3 =  -DeltaP(3)/tempDelP;
                   tempP4 =  (-DeltaP(4)-DeltaP(1)*DeltaP(4)+DeltaP(2)*DeltaP(3))/tempDelP;
                   tempP5 =  (-DeltaP(5)-DeltaP(4)*DeltaP(5)+DeltaP(3)*DeltaP(6))/tempDelP;
                   tempP6 =  (-DeltaP(6)-DeltaP(1)*DeltaP(6)+DeltaP(2)*DeltaP(5))/tempDelP;

                   tempMatrix = [1+P(1) P(3) P(5); P(2) 1+P(4) P(6); 0 0 1]*...
                       [1+tempP1 tempP3 tempP5; tempP2 1+tempP4 tempP6; 0 0 1];

                   P1 = tempMatrix(1,1)-1;
                   P2 = tempMatrix(2,1);
                   P3 = tempMatrix(1,2);
                   P4 = tempMatrix(2,2)-1;
                   P5 = tempMatrix(1,3);
                   P6 = tempMatrix(2,3);
                   P = [P1 P2 P3 P4 P5 P6]';
               end

               % Find region for g
               [tempCoordy, tempCoordx] = meshgrid(y(1):y(3),x(1):x(3));
               
               tempCoordx = tempCoordx - x0*ones(winsize+1,winsize+1);
               tempCoordy = tempCoordy - y0*ones(winsize+1,winsize+1);
               u22 = (1+P(1))*tempCoordx + P(3)*tempCoordy + (x0+P(5))*ones(winsize+1,winsize+1);
               v22 = P(2)*tempCoordx + (1+P(4))*tempCoordy + (y0+P(6))*ones(winsize+1,winsize+1);
               
               tempg = zeros(size(tempf,1)*size(tempf,2),1);
               
               [tempCoordy, tempCoordx] = meshgrid(1:winsize+1,1:winsize+1);
               tempCoordx = tempCoordx(:); tempCoordy = tempCoordy(:);
               
               parfor tempij = 1:size(tempCoordx,1)
                   tempg(tempij)= ...
                       fungInterpolation_g(u22(tempCoordx(tempij),tempCoordy(tempij)), v22(tempCoordx(tempij),tempCoordy(tempij)), ...
                       g(floor(u22(tempCoordx(tempij),tempCoordy(tempij)))-1:floor(u22(tempCoordx(tempij),tempCoordy(tempij)))+2, ...
                       floor(v22(tempCoordx(tempij),tempCoordy(tempij)))-1:floor(v22(tempCoordx(tempij),tempCoordy(tempij)))+2));
               end
               
               tempg = reshape(tempg, winsize+1, winsize+1);
               
               meang = mean(tempg(:));
               bottomg = sqrt((length(tempg(:))-1)*var(tempg(:)));
               
            otherwise
         end
        
        % ============ End of Levenberg-Marquardt method ============ 

        % Assemble b vector
        b = zeros(6,1);
        [tempCoordy, tempCoordx] = meshgrid(y(1):y(3),x(1):x(3));
        tempCoordx = tempCoordx(:); tempCoordy = tempCoordy(:);
 
        % Compute for all the pixels in the master part
        for tempij = 1:size(tempCoordx,1)
            if PassCrackOrNot == 1 % Crack pass through this local subset
                % if x0 < CrackTip(1) % Left part is master, and right part will be discarded
                if CrackPathCen(1)*x0 + CrackPathCen(2)*y0 + 1 > 0  % Bottom part is master, and top part will be discarded
                    % if tempCoordx(tempij) > CrackTip(1) && tempCoordy(tempij) < CrackTip(2)
                    %     tempCoordx(tempij) = tempCoordx(tempij)-winsize;
                    % end
                    % if tempCoordx(tempij)+1-x(1)>1 && tempCoordx(tempij)>x0-0.8*winsize
                    %     b = b+ bottomf*([DfDx(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty) DfDy(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty)]*...
                    %         [tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1 0; 0 tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1])'* ...
                    %         ((tempf(tempCoordx(tempij)+1-x(1), tempCoordy(tempij)+1-y(1))-meanf)/bottomf - ...
                    %         (tempg(tempCoordx(tempij)+1-x(1), tempCoordy(tempij)+1-y(1))-meang)/bottomg);
                    % end
                    if tempCoordx(tempij)*CrackPathCen(1) + tempCoordy(tempij)*CrackPathCen(2) + 1 < 0
                        tempCoordy(tempij) = tempCoordy(tempij)-winsize;
                    end
                    if (tempCoordy(tempij)+1-y(1) > 1) && (tempCoordy(tempij)+1-y(1) < size(tempf,2)) ...
                        && (tempCoordx(tempij)*CrackPath2(1) + tempCoordy(tempij)*CrackPath2(2) + 1 > 0)
                           b = b+ bottomf*([DfDx(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty) DfDy(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty)]*...
                            [tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1 0; 0 tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1])'* ...
                            ((tempf(tempCoordx(tempij)+1-x(1), tempCoordy(tempij)+1-y(1))-meanf)/bottomf - ...
                            (tempg(tempCoordx(tempij)+1-x(1), tempCoordy(tempij)+1-y(1))-meang)/bottomg);
                    end    
                % elseif (x0 >= CrackTip(1)) % Right part is master, and left part will be discarded
                elseif CrackPathCen(1)*x0 + CrackPathCen(2)*y0 + 1 <= 0  % Top part is master, and bottom part will be discarded
                    % if tempCoordx(tempij) < CrackTip(1) && tempCoordy(tempij) < CrackTip(2)
                    %     tempCoordx(tempij) = tempCoordx(tempij)+winsize;
                    % end
                    % if tempCoordx(tempij)<size(DfDx,1) && tempCoordx(tempij)+1-x(1)<size(tempf,1) && tempCoordx(tempij)<x0+0.5*winsize
                    %     b = b+ bottomf*([DfDx(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty) DfDy(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty)]*...
                    %         [tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1 0; 0 tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1])'* ...
                    %         ((tempf(tempCoordx(tempij)+1-x(1), tempCoordy(tempij)+1-y(1))-meanf)/bottomf - ...
                    %         (tempg(tempCoordx(tempij)+1-x(1), tempCoordy(tempij)+1-y(1))-meang)/bottomg);
                    % end
                    if tempCoordx(tempij)*CrackPathCen(1) + tempCoordy(tempij)*CrackPathCen(2) + 1 > 0
                        tempCoordy(tempij) = tempCoordy(tempij)+winsize;
                    end
                    if (tempCoordy(tempij)-DfDxStarty < size(DfDy,2)) && (tempCoordy(tempij)+1-y(1) < size(tempf,2)) ...
                        && (tempCoordx(tempij)*CrackPath1(1) + tempCoordy(tempij)*CrackPath1(2) + 1 < 0)
                            b = b + bottomf*([DfDx(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty) DfDy(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty)]*...
                            [tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1 0; 0 tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1])'* ...
                            ((tempf(tempCoordx(tempij)+1-x(1), tempCoordy(tempij)+1-y(1))-meanf)/bottomf - ...
                            (tempg(tempCoordx(tempij)+1-x(1), tempCoordy(tempij)+1-y(1))-meang)/bottomg);
                    end
                end
            else
                b = b + bottomf*([DfDx(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty) DfDy(tempCoordx(tempij)-DfDxStartx,tempCoordy(tempij)-DfDxStarty)]*...
                    [tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1 0; 0 tempCoordx(tempij)-x0 0 tempCoordy(tempij)-y0 0 1])'* ...
                    ((tempf(tempCoordx(tempij)+1-x(1), tempCoordy(tempij)+1-y(1))-meanf)/bottomf - ...
                    (tempg(tempCoordx(tempij)+1-x(1), tempCoordy(tempij)+1-y(1))-meang)/bottomg);
            end
        end
        
        tempb = b(5:6)*2/(bottomf^2)  + [mu*(P(5)-UOld(1)-vdual(1)); mu*(P(6)-UOld(2)-vdual(2))];
        %%%%%%%%%%%%%%%% Tried below, not succeed %%%%%%%%%%%%%%%%%%%
        % if CrackOrNot == 0
        %     tempb = b(5:6)*2/(bottomf^2)  + [mu*(P(5)-UOld(1)-vdual(1)); mu*(P(6)-UOld(2)-vdual(2))];
        % else
        %     tempb = b*2/(bottomf^2) +[beta*(P(1)-FOld(1)-udual(1)); beta*(P(2)-FOld(2)-udual(2));
        %                               beta*(P(3)-FOld(3)-udual(3)); beta*(P(4)-FOld(4)-udual(4));
        %                               mu*(P(5)-UOld(1)-vdual(1)); mu *(P(6)-UOld(2)-vdual(2))];
        % end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
        
        delta = 1e-2; DeltaP = [0 0 0 0 0 0];
        H2 = H(5:6,5:6)*2/(bottomf^2) + [mu 0; 0 mu];
        tempH = (H2 + delta*max(diag(H2))*eye(2));
        
        DeltaP(5:6) = -tempH\tempb;
        %%%%%%%%%%%%%%%% Tried below, not succeed %%%%%%%%%%%%%%%%%%%
        % if CrackOrNot == 0
        %     DeltaP(5:6) = -tempH\tempb;
        % else
        %     DeltaP = -(H2 + delta*diag(diag(H2))) \ tempb;
        % end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        normOfWOld = normOfWNew;
        normOfWNew = norm(tempb(:));
        
        
        if stepwithinwhile == 1
            normOfWNewInit = normOfWNew;
        end
        
        normOfWNew = normOfWNew/normOfWNewInit;
        
        if normOfWNew < tol
            break
        else
            
            tempDelP =  ((1+DeltaP(1))*(1+DeltaP(4)) - DeltaP(2)*DeltaP(3));
            if (tempDelP ~= 0)
                tempP1 =  (-DeltaP(1)-DeltaP(1)*DeltaP(4)+DeltaP(2)*DeltaP(3))/tempDelP;
                tempP2 =  -DeltaP(2)/tempDelP;
                tempP3 =  -DeltaP(3)/tempDelP;
                tempP4 =  (-DeltaP(4)-DeltaP(1)*DeltaP(4)+DeltaP(2)*DeltaP(3))/tempDelP;
                tempP5 =  (-DeltaP(5)-DeltaP(4)*DeltaP(5)+DeltaP(3)*DeltaP(6))/tempDelP;
                tempP6 =  (-DeltaP(6)-DeltaP(1)*DeltaP(6)+DeltaP(2)*DeltaP(5))/tempDelP;

                tempMatrix = [1+P(1) P(3) P(5); P(2) 1+P(4) P(6); 0 0 1]*...
                    [1+tempP1 tempP3 tempP5; tempP2 1+tempP4 tempP6; 0 0 1];

                P1 = tempMatrix(1,1)-1;
                P2 = tempMatrix(2,1);
                P3 = tempMatrix(1,2);
                P4 = tempMatrix(2,2)-1;
                P5 = tempMatrix(1,3);
                P6 = tempMatrix(2,3);
                P = [P1 P2 P3 P4 P5 P6]';
            else
                disp( 'Det(DeltaP)==0!' );
                break
            end

        end
    end  
end % end of while
 
U(1) = P(5); U(2) = P(6);
% F(1) = P(1); F(2) = P(2); F(3) = P(3); F(4) = P(4);


if (normOfWNew < tol)
else
    stepwithinwhile = 0;
end
 
if (isnan(normOfWNew)==1)
    stepwithinwhile = 0;
end

% if (elementsLocalMethodConvergeOrNot == 0)
% 	[row] = find(KappaStore == min(KappaStore(1:min(10,stepwithinwhile-1))));
% 	U(1) = PStore(row,5); U(2) = PStore(row,6);
% 	F(1) = PStore(row,1); F(2) = PStore(row,2); F(3) = PStore(row,3); F(4) = PStore(row,4);
% end

end


 