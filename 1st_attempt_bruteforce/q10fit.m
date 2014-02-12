function aff_params = q10fit()
T = [9.8 13.9 23.5 29.9 34.9 35.0 37.1]';
na_q10_factors = [.0845 .219 .566 .777 .995 1.000 1.109]';
k_q10_factors = [.0988 .165 .610 .791 .995 1.000 1.105]';

aff_params = fminsearch(@(ab) aff_errorfcn(T, na_q10_factors, ab), [1 1]);

plot(T, [na_q10_factors  k_q10_factors aff_params(1)*T+aff_params(2)])
legend('Na^+', 'K^+', 'Na^+ linear fit', 'Location', 'Northwest')


function error = aff_errorfcn(pts, ypts, params)
pred = params(1)*pts + params(2);
sqerr = (ypts - pred).^2;
error = sum(sqerr);

