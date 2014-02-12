%cytoplasmic resistivity R_i

T = [7.7 9.8 13.0 23.5 29.9 34.9 35.0 37.1]';

q10_i = .8;
ref_R_i = 140; %Ohms*cm
ref_temp_i = 36; %C

calc_R_i = ref_R_i./(q10_i.^((ref_temp_i-T)/10));
calc_R_i = round(calc_R_i*10)/10;

given_R_i = [263.3 251.2 234.9 185.0 160.4 143.5 143.2 136.6]';

mean_sq_error = mean((calc_R_i-given_R_i).^2);

%leak conductance
ref_res = 10; %Kohm/cm^2
ref_g_l = 1/mem_res;
ref_temp_r = 35; 

q10_r = 1.85;

res_fac = ref_res./(q10_r.^((ref_temp_r-T)/10));

g_l_fac = 1./res_fac;

fprintf(['Temperature     Given R_i (Fohlmeister)    Calculated R_i    ' ...
         'Calculated conductance (mS/cm^2)\n'])

for n = numel(T):-1:1
fprintf('  %4.1f                %4.1f                   %4.1f                      %4.5f\n',...
          T(n),given_R_i(n),calc_R_i(n), g_l_fac(n));
end