function [mod] = compfunc_ultim(RATES)

m = 0;

global t50 A0 G0 R0 K0

A0 = RATES(1); 
t50 = RATES(2);
G0 = RATES(3);
K0 = RATES(4);

x = linspace(0,max_val,1000);
t = linspace(0,1,20);

if K0 >= 0

options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
sol = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc,x,t, options);


% sol = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc,x,t);

% Extract the first solution component as mod.
mod = sol(20,:,1);

else
mod = NaN(1,1000);
end

%--------------------------------------------------
% --------------------------------------------------------------
function [c,f,s] = pdex1pde(x,t,u,DuDx)
global t50 A0 G0 K0

u = max(u, 0);
nu  = A0./(1+exp(-G0*(x-t50)));
c = 1;
f =  nu*u; 
s = -3*nu*u/x - K0*u; % 
% s = - K0*u; % 

% --------------------------------------------------------------
function u0 = pdex1ic(x)
u0 = init_pdf_N2(x);
% Add small regularization to prevent exactly zero values
u0 = u0 + 1e-10;

% --------------------------------------------------------------
function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t)
global t50 A0 G0

pl = ul;
ql = 0;

pr = ur;
qr = 1e-10;

