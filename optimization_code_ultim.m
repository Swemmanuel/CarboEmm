clc
close all
clear all
 
% Code finds the optimal rate constants to fit dataset

% Initial guess for size dependent rate
% for dissolution, A0 is -ve
% for precipitation, A0 is +ve
% y = A0./(1+exp(-G0*(x-t50))); 20.50	67.29	0.08	0.73

% Initial guess for required duration of simulation
         
A_0 = 1.3699; 
T50 = 31.1404;
G_0 = -0.0251 ; 
K_0 = 0.4407;


RATE = [A_0, T50, G_0, K_0];

% evaluate model for first guess  
[model] = compfunc_ultim(RATE);

% call datafile for visualization
xD = linspace(0, max_val*2, 1000); % factor of 2 is to convert radius into diameter
xr = xD/2; 
y0 = init_pdf_N2(xr);
yf = fin_pdf_N2(xr);

total_area = trapz(xr, y0);
fprintf('Total porosity, initial curve: %.4f\n', total_area);

total_area = trapz(xr, yf);
fprintf('Total porosity, final curve: %.4f\n', total_area);

% plot data
figure(1)
firstcolor = [31/255, 119/255, 180/255];  % [0.8500 0.3250 0.0980];
plot(xD, y0/2,'-', 'LineWidth', 2, 'Color', firstcolor)
hold on
secondcolor = [255/255, 127/255, 14/255]; 
plot(xD, yf/2, 'LineWidth', 2.5, 'Color', secondcolor)
axis([0, max_val*2, 0, 0.0052])
xlabel('pore diameter [\mum]', 'fontsize', 14)
ylabel('\phi_r [\mum^{-1}]', 'fontsize', 14)

legend('boxoff') 
box on

% calculate error
initial_error = diffnorm_pores_ultim(RATE)

% optimization code
OPTS = optimset('TolFun', 1e-12, 'Tolx', 1e-12, 'MaxIter', 100, 'display', 'iter');
[new_RATE err ] = fminsearch('diffnorm_pores_ultim', RATE, OPTS)

%  optimized rates!
RATE = new_RATE;

% % calculate optimized model  
[optim_model] = compfunc_ultim(RATE);
    result = optim_model;
    result(result < 0) = 0;
    optim_model = result;


NEW_RATE = RATE; 

A0 = NEW_RATE(1);
t50 = NEW_RATE(2);
G0 = NEW_RATE(3);
K0 = NEW_RATE(4)

% plot optimized model
figure(1)
hold on
ccolor = [0.6 0.6 0.6];
plot(xD, optim_model/2, 'LineWidth', 2, 'Color', ccolor)
box on
h_legend = legend('before cementation', 'after cementation', 'optimized model');
set(h_legend,'FontSize',14)



