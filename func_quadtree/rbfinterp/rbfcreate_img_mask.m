function options = rbfcreate_img_mask(x, y, ImgMask, varargin)
%RBFCREATE Creates an RBF interpolation
%   OPTIONS = RBFSET(X, Y, 'NAME1',VALUE1,'NAME2',VALUE2,...) creates an   
%   radial base function interpolation 
%   
%   RBFCREATE with no input arguments displays all property names and their
%   possible values.
%   
%RBFCREATE PROPERTIES
% 

%
% Alex Chirokov, alex.chirokov@gmail.com
% 16 Feb 2006
tic;
% Print out possible values of properties.
if (nargin == 0) & (nargout == 0)
  fprintf('               x: [ dim by n matrix of coordinates for the nodes ]\n');
  fprintf('               y: [   1 by n vector of values at nodes ]\n');
  fprintf('     RBFFunction: [ gaussian  | thinplate | cubic | multiquadrics | {linear} ]\n');
  fprintf('     RBFConstant: [ positive scalar     ]\n');
  fprintf('       RBFSmooth: [ positive scalar {0} ]\n');
  fprintf('           Stats: [ on | {off} ]\n');
  fprintf('\n');
  return;
end
Names = [
    'RBFFunction      '
    'RBFConstant      '
    'RBFSmooth        '
    'Stats            '
];
[m,n] = size(Names);
names = lower(Names);

options = [];
for j = 1:m
  options.(deblank(Names(j,:))) = [];
end

%**************************************************************************
%Check input arrays 
%**************************************************************************
[nXDim nXCount]=size(x);
[nYDim nYCount]=size(y);

if (nXCount~=nYCount)
  error(sprintf('x and y should have the same number of rows'));
end;

if (nYDim~=1)
  error(sprintf('y should be n by 1 vector'));
end;

options.('x')           = x;
options.('y')           = y;
%**************************************************************************
%Default values 
%**************************************************************************
options.('RBFFunction') = 'linear';
options.('RBFConstant') = (prod(max(x')-min(x'))/nXCount)^(1/nXDim); %approx. average distance between the nodes 
options.('RBFSmooth')   = 0;
options.('Stats')       = 'off';

%**************************************************************************
% Argument parsing code: similar to ODESET.m
%**************************************************************************

i = 1;
% A finite state machine to parse name-value pairs.
if rem(nargin-3,2) ~= 0
  error('Arguments must occur in name-value pairs.');
end
expectval = 0;                          % start expecting a name, not a value
while i <= nargin-3
  arg = varargin{i};
    
  if ~expectval
    if ~isstr(arg)
      error(sprintf('Expected argument %d to be a string property name.', i));
    end
    
    lowArg = lower(arg);
    j = strmatch(lowArg,names);
    if isempty(j)                       % if no matches
      error(sprintf('Unrecognized property name ''%s''.', arg));
    elseif length(j) > 1                % if more than one match
      % Check for any exact matches (in case any names are subsets of others)
      k = strmatch(lowArg,names,'exact');
      if length(k) == 1
        j = k;
      else
        msg = sprintf('Ambiguous property name ''%s'' ', arg);
        msg = [msg '(' deblank(Names(j(1),:))];
        for k = j(2:length(j))'
          msg = [msg ', ' deblank(Names(k,:))];
        end
        msg = sprintf('%s).', msg);
        error(msg);
      end
    end
    expectval = 1;                      % we expect a value next
    
  else
    options.(deblank(Names(j,:))) = arg;
    expectval = 0;      
  end
  i = i + 1;
end

if expectval
  error(sprintf('Expected value for property ''%s''.', arg));
end

    
%**************************************************************************
% Creating RBF Interpolatin
%**************************************************************************

switch lower(options.('RBFFunction'))
      case 'linear'          
        options.('rbfphi')   = @rbfphi_linear;
      case 'cubic'
        options.('rbfphi')   = @rbfphi_cubic;
      case 'multiquadric'
        options.('rbfphi')   = @rbfphi_multiquadrics;
      case 'thinplate'
        options.('rbfphi')   = @rbfphi_thinplate;
      case 'gaussian'
        options.('rbfphi')   = @rbfphi_gaussian;
    otherwise
        options.('rbfphi')   = @rbfphi_linear;
end

phi       = options.('rbfphi');

A=rbfAssemble(x, phi, ImgMask, options.('RBFConstant'), options.('RBFSmooth'));

b=[y'; zeros(nXDim+1, 1)];                       

%inverse
rbfcoeff=A\b;

%SVD
% [U,S,V] = svd(A);
% 
% for i=1:1:nXCount+1
%     if (S(i,i)>0) S(i,i)=1/S(i,i); end;   
% end;    
% rbfcoeff = V*S'*U*b;


options.('rbfcoeff') = rbfcoeff;


if (strcmp(options.('Stats'),'on'))
    fprintf('%d point RBF interpolation was created in %e sec\n', length(y), toc);  
    fprintf('\n');
end;

function [A]=rbfAssemble(x, phi, ImgMask, const, smooth)
[dim n]=size(x);
A=zeros(n,n);
for i=1:n
    for j=1:i

        %%%%% Whether i and j are connected %%%%%
%         Img_IJ = func_DrawLine(zeros(size(ImgMask)), round(x(1,i)), round(x(2,i)), round(x(1,j)), round(x(2,j)), 1);
%         Img_IJ_connected_or_not = sum(sum(Img_IJ .* (1-double(ImgMask))));
        
        Img_IJ_Path_x = floor(min(x(1,[i,j]))) : 1 : ceil(max(x(1,[i,j])));
        if length(Img_IJ_Path_x) > 1
            Img_IJ_Path_y = round(interp1([x(1,i),x(1,j)],[x(2,i),x(2,j)],Img_IJ_Path_x));
        else
            Img_IJ_Path_y = floor(min(x(2,[i,j]))) : 1 : ceil(max(x(2,[i,j])));
            if length(Img_IJ_Path_y) > 1
                Img_IJ_Path_x = round(interp1([x(2,i),x(2,j)],[x(1,i),x(1,j)],Img_IJ_Path_y));
            end
        end
        % figure, plot(Img_IJ_Path_x,Img_IJ_Path_y,'.');
        indImgMask = sub2ind(size(ImgMask),Img_IJ_Path_x,Img_IJ_Path_y);
        Img_IJ_connected_or_not = sum( 1-ImgMask(indImgMask) );
        
        if Img_IJ_connected_or_not < 1
        
            r=norm(x(:,i)-x(:,j));
            temp=feval(phi,r, const); 
            A(i,j)=temp;
            A(j,i)=temp;
            
        end
         
    end
    A(i,i) = A(i,i) - smooth;
end
% Polynomial part
P = [ones(n,1) x'];
A = [ A      P
      P' zeros(dim+1,dim+1)];

%**************************************************************************
% Radial Base Functions
%************************************************************************** 
function u=rbfphi_linear(r, const)
u=r;

function u=rbfphi_cubic(r, const)
u=r.*r.*r;

function u=rbfphi_gaussian(r, const)
u=exp(-0.5*r.*r/(const*const));

function u=rbfphi_multiquadrics(r, const)
u=sqrt(1+r.*r/(const*const));

function u=rbfphi_thinplate(r, const)
u=r.*r.*log(r+1);