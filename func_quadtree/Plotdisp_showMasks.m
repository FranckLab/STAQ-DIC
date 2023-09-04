function Plotdisp_showMasks(U,coordinatesFEMWorld,elementsFEM,CurrentImgMask,varargin)
%PLOTDISP_SHOW: to plot DIC solved displacement components
%   Plotdisp_show(U,coordinatesFEM,elementsFEM)
% ----------------------------------------------
%
%   INPUT: U                 Displacement vector: 
%                            U = [Ux_node1, Uy_node1, Ux_node2, Uy_node2, ... , Ux_nodeN, Uy_nodeN]';
%          coordinatesFEM    FE mesh coordinates
%          elementsFEM       FE mesh elements
%          DICpara           chosen DIC parameters
%          EdgeColorOrNot    show edge color or not
%
%   OUTPUT: Plots of x-displacement field and y-displacement field.
%
%   TODO: users could change caxis range based on their own choices.
%
% ----------------------------------------------
% Author: Jin Yang.  
% Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
% Last date modified: 2020.12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Initialization
warning off;  U = full(U);
%run('./img_Vito/Black_rainbow.m');

%%%%% Parse Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[DICpara,EdgeColorOrNot] = parseargs(varargin);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% convert pixel unit to the physical world unit %%%%%
try um2px = DICpara.um2px; 
catch um2px = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try EdgeColorOrNot = EdgeColorOrNot;
catch EdgeColorOrNot = 'EdgeColor';
end

Image2PlotResults = 0;
disp_u = U(1:2:end); disp_v = U(2:2:end);
coordinatesFEMWorldDef = [coordinatesFEMWorld(:,1)+Image2PlotResults*disp_u, coordinatesFEMWorld(:,2)+Image2PlotResults*disp_v];

%%%%%%%%%%% JY!!!Mask START %%%%%%%%%%%%%%%
if ~isempty(CurrentImgMask)
    for tempi = 1:size(coordinatesFEMWorldDef,1)
%         try
            if CurrentImgMask( floor(coordinatesFEMWorldDef(tempi,1)/um2px), ...
                    (size(CurrentImgMask,2)+1-ceil(coordinatesFEMWorldDef(tempi,2)/um2px)) ) == 0
                coordinatesFEMWorldDef(tempi,:) = [nan,nan];
            end
%         catch
%             coordinatesFEMWorldDef(tempi,:) = [nan,nan];
%         end     
    end
end
%%%%%%%%%%% JY!!!Mask END %%%%%%%%%%%%%%%


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 1) dispx u ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; show([],elementsFEM(:,1:4),coordinatesFEMWorld,U(1:2:end),EdgeColorOrNot); 

title('$x-$displacement $u$','FontWeight','Normal','Interpreter','latex');
view(2); set(gca,'fontsize',18); axis tight; axis equal; colorbar; % view([90 -90])
if um2px==1, xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
else, xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex');
end
set(gcf,'color','w'); colormap jet; box on;
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';

%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet);  % D Sample 
% caxis([-0.1,0.1]) % foam
% caxis([-0.004,0.004]); % Sample 12
% caxis([0,0.1]); 
% colormap(black_rainbow_plus); caxis([-0.56,0.56]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====== 2) dispx v ======
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; show([],elementsFEM(:,1:4),coordinatesFEMWorld,U(2:2:end),EdgeColorOrNot); 

title('$y-$displacement $v$','FontWeight','Normal','Interpreter','latex');
view(2); set(gca,'fontsize',18); axis tight; axis equal; colorbar;
if um2px==1, xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
else, xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex');
end
set(gcf,'color','w'); colormap jet; box on;
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';

%%%%%% TODO: manually modify colormap and caxis %%%%%%
% colormap(jet);  % D Sample 
% caxis([-0.1,0.1]) % foam
% caxis([-0.004,0.004]); % Sample 12
% caxis([0,0.1]); 
% colormap(black_rainbow_plus); caxis([-0.12,0.12]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function [DICpara,EdgeColorOrNot] = parseargs(vargin)
DICpara=[]; EdgeColorOrNot=[];  
 
try 
    DICpara=vargin{1};
    try
        EdgeColorOrNot=vargin{2};
    catch
    end
catch
    
end

end
