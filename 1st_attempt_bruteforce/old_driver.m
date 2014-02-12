%% driver for temperature comparison in bf4

g_l_fcn = @(an_x_grid) (an_x_grid)*2.5 + .5;
l = .1;

Ntemp = [1 2 5 10 20];
Nbf = 10;
bf_mod_no = 4;

x_vals = linspace(0,l,Nbf)';
true_g_l = g_l_fcn(x_vals);

xstep = zeros(Nbf*2+1, 1);
out_xvals = linspace(0,l,Nbf+1)';
xstep(1:2:2*Nbf) = out_xvals(1:end-1)';
xstep(2:2:2*Nbf) = out_xvals(2:end)'; 

step_g_l = zeros(Nbf*2, 1);

if 0

NNt = numel(Ntemp);
NTr = numel(Trange);

error = zeros(NNt, NTr);
g_l_rec = zeros(Nbf, NTr, NNt);



for n = 1:NNt

      fprintf('\nNtemp = %d, Trange = %d\n', Ntemp(n), Trange(m));
      figure(1)
      [g_l_rec(:,m,n) error(n,m)] = bf4('step', Ntemp(n), ...
                                        Trange(m), Nbf, bf_mod_no);

      clear persistent
      figure(6)
      for b = 1:Nbf
          step_g_l(n*2-1:n*2) = g_l_rec(b,m,n);
      end
      subplot(NTr, NNt, (m-1)*NNt + n)
      plot(x_vals, true_g_l, xstep(1:Nbf*2), step_g_l);
      title(['Ntemp = ' num2str(Ntemp(n)) ', Range = ' ...
              num2str(Trange(m)) ]);


end


end

figure(2)


for n = 1:NNt
        for b = 1:Nbf
            step_g_l(b*2-1:b*2) = g_l_rec(b,m,n);
        end
        subplot(1, NNt, n)
        plot(x_vals, g_l_fcn(x_vals), xstep(1:Nbf*2), step_g_l);
        title(['Ntemp = ' num2str(Ntemp(n)) '])
end


figure(3)
subplot 121
    plot(error)
    xlabel('Temperature Range', 'FontSize', 15)

    surf(error);
    ylabel('Number of temperatures', 'FontSize', 15)
    xlabel('Temperature Range', 'FontSize', 15)
end
title('Error in g_l: semicontinuous', 'FontSize', 18);


%alternate error measure: just centers of basis functions
l = .1;

rep_true_g_l = repmat(true_g_l, [1, numel(Trange), numel(Ntemp)]);
sum_sq_error = sum((g_l_rec - rep_true_g_l).^2);
plottable_node_error = reshape(sum_sq_error, numel(Trange),numel(Ntemp))';


subplot 122
if numel(Trange)<2
    plot(plottable_node_error)
    xlabel('Number of temperatures', 'FontSize', 15)
elseif numel(Ntemp)<2
    plot(plottable_node_error)
    xlabel('Temperature Range', 'FontSize', 15)
else
    surf(plottable_node_error);
    ylabel('Number of temperatures', 'FontSize', 15)
    xlabel('Temperature Range', 'FontSize', 15)
end

title('Error in g_l: at basis nodes', 'FontSize', 18);



