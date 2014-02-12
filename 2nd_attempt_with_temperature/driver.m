%driver

global cond_vec
Nexp = 10;
Nx = 140;
cond_vec = zeros(Nx, Nexp,2);

for n = 1:Nexp
   fifth_tempx(n) 
end

save cond_vec cond_vec

max_y = 1e13;

figure(99)
 subplot 121
cmap = colormap('jet');
inc = floor(size(jet,1)/Nexp);
cvh = semilogy(cond_vec(:,:,1)); axis([0 140 0 max_y])
for n = 1:Nexp
    set(cvh(n), 'Color', cmap(n*inc,:));
end
title('Condition number of estimate A', 'FontSize' , 16)
xlabel('Iterations (space step)')

 subplot 122
known_cvh = semilogy(cond_vec(:,:,2)); axis([0 140 0 max_y])

for n = 1:Nexp
    set(known_cvh(n), 'Color', cmap(n*inc,:));
end
title('Condition number of true A', 'FontSize' , 16)
xlabel('Iterations (space step)')
legend('5', '10', '15', '20', '25', '30', '35', '40', '45', '50')

figure(100)
subplot 121
cvh3 = semilogy(cond_vec(:,:,3)); axis([0 140 0 max_y])
for n = 1:Nexp
    set(cvh3(n), 'Color', cmap(n*inc,:));
end
title('cond. number of estimate checkmat')
subplot 122
cvh4 = semilogy(cond_vec(:,:,4)); axis([0 140 0 max_y])
for n = 1:Nexp
    set(cvh4(n), 'Color', cmap(n*inc,:));
end
title('cond. no of true checkmat')
legend('5', '10', '15', '20', '25', '30', '35', '40', '45', '50')

