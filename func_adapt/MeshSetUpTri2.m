%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: FEM Kuhn Simplex Triangle Mesh Set Up
% Author: Jin Yang
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [coordinatesFEM,elementsFEM,coordinates,elements,dirichlet,neumann,x,y,xCen,yCen,M,N] = MeshSetUpTri2(x,y,xCen,yCen,winstepsize)

% ========= mesh for global method ===========
M = size(x,2);  N = size(x,1); % N is vertically in image; M is horizontally in image;
coordinatesFEM = zeros(M*N,2);

x = x'; y = y';
% I have transpose x and y because Matlab is read matrix in column direction
for tempi = 1:(M*N)
    coordinatesFEM(tempi,:) = [x(tempi),y(tempi)];
    % x is horizontal position in the image
    % y is vertical position in the image
end


elementsFEM = zeros(2*(M-1)*(N-1),3); % The factor 4 is because I cut one Q4 element into 4 Kuhn triangle elements
% Triangle elements have the following structures.
%    3-----4 
%    |    /|        
%    |   / |        
%    |  /  |
%    | /   |
%    |/    |
%    1-----2
% {4-1-2}; {1-4-3};  

for j = 1:N-1
    for i = 1:M-1
        
        Node1=(j-1)*(M)+i;  Node2=(j-1)*(M)+i+1;    
        Node3=j*(M)+i;      Node4=j*(M)+i+1;     
        elementsFEM(2*((j-1)*(M-1)+i)-1,:) = [Node4, Node1, Node2];
        elementsFEM(2*((j-1)*(M-1)+i)-0,:) = [Node1, Node4, Node3];
  
    end
end


% ========== mesh for local method ===========
% N is vertically in image; M is horizontally in image;
coordinates = zeros((M+1)*(N+1)+M*N,2);
gridxtemp = [x(1,1)-0.5*winstepsize; x(:,1)+0.5*winstepsize];
gridytemp = [y(1,1)-0.5*winstepsize y(1,:)+0.5*winstepsize]';
clear gridx gridy
[gridx,gridy] = meshgrid(gridxtemp,gridytemp);
gridx = gridx'; gridy = gridy';

for tempi = 1:((M+1)*(N+1))
    coordinates(tempi,:)  = [gridx(tempi),gridy(tempi)];
end
for tempi = 1:(M*N)
    coordinates((M+1)*(N+1)+tempi,:) = [x(tempi),y(tempi)];
end

elements = zeros( M*N + ((M-1)*(N-1)) ,4);
for j = 1:N
    for i = 1:M
        elements((j-1)*(M)+i, :) = [(j-1)*(M+1)+i, (j-1)*(M+1)+i+1, j*(M+1)+i+1, j*(M+1)+i];
    end
end
for j = 1:N-1
    for i = 1:M-1
        elements(M*N+(j-1)*(M-1)+i, :) = [(M+1)*(N+1)+(j-1)*(M)+i, (M+1)*(N+1)+(j-1)*(M)+i+1, ...
                                            (M+1)*(N+1)+j*(M)+i+1, (M+1)*(N+1)+j*(M)+i];
    end
end

 
%% % ------ Assign Boundary conditions ------
% neumann = [1:1:M, 1:M:M*N, M:M:M*N, (M*(N-1)+1):1:M*N]';
% neumann = unique(neumann);
% dirichlet = [];

% ======== Assign BC values ==========
% -------- dirichlet BC --------
dirichlet = [linspace(1,(M-1),(M-1))', linspace(2,M,(M-1))' ;
            linspace(1,1+M*(N-2),(N-1))', linspace((1+M),(1+M*(N-1)),(N-1))' ;
            linspace(M,M*(N-1),(N-1))', linspace(2*M,M*N,(N-1))' ;
            linspace((M*(N-1)+1), (M*N-1), (M-1))', linspace((M*(N-1)+2), (M*N), (M-1))' ];
%dirichlet = [];
% -------- neumann BC --------       
neumann = [linspace(1,(M-1),(M-1))', linspace(2,M,(M-1))', zeros(M-1,1), -ones(M-1,1) ;
             linspace(1,1+M*(N-2),(N-1))', linspace((1+M),(1+M*(N-1)),(N-1))', -ones(N-1,1), zeros(N-1,1) ;
             linspace(M,M*(N-1),(N-1))', linspace(2*M,M*N,(N-1))', ones(N-1,1), zeros(N-1,1) ;
             linspace((M*(N-1)+1), (M*N-1), (M-1))', linspace((M*(N-1)+2), (M*N), (M-1))', zeros(M-1,1), ones(M-1,1) ];
%neumann = [];


