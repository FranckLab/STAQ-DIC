function Plotdisp_showQuadtree(U,coordinatesFEM,elementsFEM)
%FUNCTION Plotdisp_showQuadtree(U,coordinatesFEM,elementsFEM)
% To plot DIC solved displacement components
% ----------------------------------------------
%
%   INPUT: U                 Displacement vector: U = [Ux_node1, Uy_node1, Ux_node2, Uy_node2, ... , Ux_nodeN, Uy_nodeN]';
%          coordinatesFEM    FE mesh coordinates
%          elementsFEM       FE mesh elements
%
%   OUTPUT: Plots of x-displacement field and y-displacement field.
%
% ----------------------------------------------
% Author: Jin Yang.  
% Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
% Last time updated: 11/2020.
% ==============================================


%% Initialization
warning off;  U = full(U);

%% ====== 1) dispx u ======
figure; % show([],elementsFEM(:,1:4),coordinatesFEM,U(1:2:end) ); 
showQuadtree(elementsFEM(:,1:4),coordinatesFEM,U(1:2:end));

title('$x-$displacement $u$','FontWeight','Normal','Interpreter','latex');
view(2); set(gca,'fontsize',18); axis tight; axis equal; colorbar; % view([90 -90])

xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
set(gcf,'color','w'); colormap jet;
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';

colormap jet; box on;


%% ====== 2) dispy v ======
figure; % show([],elementsFEM(:,1:4),coordinatesFEM,U(2:2:end) ); 
showQuadtree(elementsFEM(:,1:4),coordinatesFEM,U(2:2:end));
title('$y-$displacement $v$','FontWeight','Normal','Interpreter','latex');
view(2); set(gca,'fontsize',18); axis tight; axis equal; colorbar;

xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
set(gcf,'color','w'); colormap jet;
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';

colormap jet; box on;



