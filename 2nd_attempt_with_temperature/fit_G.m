function aff_params = fit_G()
%Ryan George
%VIGRE 2010, Rice University
% fit_G fits conductance G to data from the paper Fohlmeister et. al to an
% affine plot

%temperatures tested
T = [7.7 9.8 13.9 16.7 23.5 29.9 34.9 37.1]; %C

%axial resistivity
R_i = [263.3 251.2 234.9 185.0 160.4 143.5 143.2 136.6]/1000; %KOhm*cm

%conductivity
G = 1./R_i;

aff_params = fminsearch(@(ab) aff_errorfcn(T, G, ab), [1 1]);

%plot(T,G, T, fin_params(1)*T + fin_params(2))

function error = aff_errorfcn(pts, ypts, params)
pred = params(1)*pts + params(2);
sqerr = (ypts - pred).^2;
error = sum(sqerr);