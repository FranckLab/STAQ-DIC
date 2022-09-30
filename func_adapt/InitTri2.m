%% Correct and Smooth initial displacements
function U0 =InitTri2(u,v,Phi,uCen,vCen,PhiCen,M,N)


% ========= Delete bad intial points ==========
% In inpaint_nans part, I use Spring Model for u and v initial guesses
uInit = u; vInit = v; 
threshod = 0.3;
[row, col] = find(Phi<(1-threshod));
uInit(row,col) = NaN; vInit(row,col) = NaN;

uInit = inpaint_nans(uInit,4);
vInit = inpaint_nans(vInit,4);
 

% ========= Set initial values for U =========
U0 = zeros(2*(M*N ),1);
uInit = uInit'; vInit = vInit'; % Here transpose bcz following give valus in column order
uCen = uCen';  vCen = vCen';
for tempi = 1:M*N  
    U0(2*tempi-1) = uInit(tempi);  
    U0(2*tempi)   = vInit(tempi);  
end
 

disp('Finish setting up mesh and assigning initial value!')
  