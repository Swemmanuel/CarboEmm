function  [err] = diffnorm_pores_ultim(RATE)
% This code calculates the error between the model and the data

% evaluate model
model = compfunc_ultim(RATE);

% call data
x = linspace(0, max_val, 1000);
datay = fin_pdf_N2(x);

if RATE(1) < 0
    err = NaN;
else
% calculate error
err = norm(datay - model);
end




