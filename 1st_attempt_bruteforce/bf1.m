function [error g_l_output] = bf1(basis_fcn, Ntemp, Trange, Ngl_els)
% Ryan George
% VIGRE, Rice University
% This code naively recovers time-integrated V from initial and Neumann
% data, then seeks to replicate g_l(x) using this data
%


a = 1; %micro_m
a = a/(1e4); %cm
C_m = 1; %microF/cm^2

l = .1; %cm
dx = 4e-3; %cm
x_grid = (0:dx:l)'; %'
Nx = numel(x_grid);

t_final = 10; %ms
dt = 3e-2; %cm from Fohlmeister: times digitized at 10, 20, 40 kHz 
  %                                            dt = .1, .05, .025 ms
time = 0:dt:t_final;
Nt = numel(time);

if nargin<2
   Ntemp = 10;
end

maxtemp = 32;
mintemp = 27;

if nargin>1
    mintemp = maxtemp - Trange;
end


Tgrid = linspace(mintemp,maxtemp, Ntemp);
%dtemp = Tgrid(2)-Tgrid(1);

g_0 = .5; %mS/cm^2
g_l_fcn = @(an_x_grid) an_x_grid*0 + g_0; %(an_x_grid .^2)*1e3 + g_0;
g_l = g_l_fcn(x_grid); %must be a vector of size(x_grid)

%for Na:
%assume G linear
aff_params = fit_G(); %KOhms
G = @(T) aff_params(1)*T +aff_params(2); %KOhms
Gprime = @(T) aff_params(1);

Gvec = G(Tgrid)*pi*(a^2)/dx;

C = 2*pi*a*dx*C_m;
tau = C_m./g_l;

%stim = zeros(Nt+1,1); %current fed into soma end of cell
stim_amp = 1e-4; %microV
mfac = (dt/2)/(2*pi*a*dx)/C_m;%(dt/2)/(2*pi*a*dx)/C_m;
% t1 = 1;
% t2 = 3;
% stim_indices = round(t1/dt): round(t2/dt);
% stim(stim_indices) = stim_amp*mfac;

% stim = (exp(-((time-2).^6)/.01)-exp(-((time-4).^6)/.01))*...
%             stim_amp*mfac;        
stim = (time>1).*(time<1.5)-(time>3).*(time<3.5)+2*(time>5).*(time<5.5)-2*(time>7).*(time<7.5);
stim(end+1)=0; %to avoid an error in the last loop
stim = stim*stim_amp*mfac;



v_tru = cab_eq(x_grid, dx, g_l, Nt, dt, G(Tgrid), stim, a, tau);

clf;
title('Forward-Constructed V', 'FontSize', 15)
mod_num = 30;
plot_t_trunc_vec = 1:mod_num:Nt;
[tg xg] = meshgrid(time(plot_t_trunc_vec), x_grid);
%plot3(tg,xg,v_tru(:,plot_t_trunc_vec,1), 'b')




%test and recover g_l
if nargin < 3
    Ngl_els = 8; %number of steps in g_l step function
end

xgl_pts = linspace(0,l,Ngl_els+1);
xgrid_to_gl_vector = zeros(Nx,1);

if nargin < 1
    basis_fcn = 'hat';
end


switch basis_fcn

    case 'step'
        for n = 1:Ngl_els
            xgrid_to_gl_vector = xgrid_to_gl_vector + (x_grid>=xgl_pts(n));
        end

        %     case 'hat'
        %         intvl = round(Nx/2/Ngl_els);
        %         for n = 0:2*Ngl_els-1
        %             xgrid_to_gl_vector(intvl*n+1:intvl*(n+1)) = n+1;
        %         end
        %         xgrid_to_gl_vector = xgrid_to_gl_vector(1:Nx);

    case 'polyn'
        xgrid_to_g1_vector = x_grid;

    case 'sigmoid'
        xgrid_to_g1_vector = x_grid;
        
end




runv = @(test_g_l) cab_eq(x_grid, dx, test_g_l, Nt, dt, G(Tgrid), stim, a, tau);

minfcn = @(g_l_info) g_l_minimization(g_l_info, x_grid,  xgrid_to_gl_vector,...
    v_tru, runv, g_l, basis_fcn);

options = optimset('MaxFunEvals', 1e3, 'TolFun', 1e-10);

tic
figure(2)
switch basis_fcn

    case 'step'
        initial_guess = ones(Ngl_els,1)*g_0/2;
        g_l_rec = fmincon(@(g_l_info) minfcn(g_l_info), initial_guess,[],[],[],[],...
            ones(Ngl_els,1)*.05,ones(Ngl_els,1)*3*max(g_l),[], options);
        initial_guess = ones(Ngl_els,1)*g_0+.05;
    case 'hat'
        initial_guess = ones(Ngl_els,1)*g_0+.05;
        %       initial_guess = g_l_fcn(linspace(0,l,Ngl_els))*1.35;
 
        g_l_rec = fmincon(@(g_l_info) minfcn(g_l_info), initial_guess,[],[],[],[],...
            ones(Ngl_els,1)*.05,ones(Ngl_els,1)*3*max(g_l),[], options);

    case 'polyn'
        initial_guess = [zeros(Ngl_els-1,1); 1];
        initial_guess = g_l(linspace(0,l,Ngl_els))+.1;

        g_l_rec = fminunc(@(g_l_info) minfcn(g_l_info), initial_guess);

    case 'sigmoid'
        initial_guess = [1 1 1 1]; %[a b h s]
        g_l_rec = fminunc(@(g_l_info) minfcn(g_l_info), initial_guess);
end

elapsed_time = toc

figure(4)
%center_pts = (xgl_pts(2:end) + xgl_pts(1:end-1))/2;


%plot

switch basis_fcn
    case 'step'
        xstep = zeros(Ngl_els*2, 1);
        step_g_l = zeros(Ngl_els*2, 1);
        for n = 1:Ngl_els
            xstep(n*2:n*2+1) = xgl_pts(n+1);
            step_g_l(n*2-1:n*2) = g_l_rec(n);
        end
        plot(x_grid, g_l, xstep(1:Ngl_els*2), step_g_l, 'g');

    case 'hat'


        intvl = round(Nx/(Ngl_els-1))*2;
        g_l_fin = zeros((Ngl_els+1)*intvl, 1);
        peak = (intvl+1)/2;
        half_width = (intvl-1)/2;
        hat_vals = 1 - abs(peak - (1:intvl)')/half_width;

        %first hat (half; falls off left boundary)
        g_l_fin(1:intvl/2) = ...
            g_l_fin(1:intvl/2) + g_l_rec(1)*hat_vals(intvl/2+1:intvl);

        %add on all other hats
        for n = 0:Ngl_els-2
            g_l_fin(intvl*n/2+1:intvl/2*(n+2)) = ...
                g_l_fin(intvl*n/2+1:intvl/2*(n+2)) + g_l_rec(n+2)*hat_vals;
        end
        g_l_fin = g_l_fin(1:Nx);
        plot(x_grid, [g_l(:) g_l_fin(:)]);

    case 'polyn'
        g_l_fin = polyval(g_l_rec, x_grid/l);
        plot(x_grid, [g_l(:) g_l_fin(:)]);
        
    case 'sigmoid'
        a = g_l_rec(1); b = g_l_rec(2); h = g_l_rec(3); s = g_l_rec(4);
        g_l_fin = b + (a-b)./(1+exp((h-x_grid*100/l)/s));
        plot(x_grid, [g_l(:) g_l_fin(:)]);
        
end

title({'Recovering Conductance: Brute Force on Lateral Overdetermination';...
    ['dx = ' num2str(dx) ', dt = ' num2str(dt) ', #' basis_fcn ...
    's = ' num2str(Ngl_els)];...
    ['\theta range = [' num2str(Tgrid(1)) ' ' num2str(Tgrid(end)) ']'...
    ', #\theta = ' num2str(Ntemp)]}, 'FontSize', 15)
axis([0 l 0 2*max(g_l)])
xlabel('Length (cm)', 'FontSize', 12)
ylabel('g_L (mS/cm^2)', 'FontSize', 12)
legend('true g_L', 'recovered g_L');


error = sum((g_l_fin(:) - g_l(:)).^2);
g_l_output = g_l_rec(:);

end



function err = g_l_minimization(g_l_info, x_grid, xgrid_to_gl_vector, v_tru, ...
    fhandle, g_l_tru, basis_fcn)

switch basis_fcn

    case 'step'
        g_l_test = zeros(size(xgrid_to_gl_vector));
        for n = 1:numel(g_l_info)
            g_l_test(xgrid_to_gl_vector==n) = g_l_info(n);
        end

        if sum(g_l_info<=0)
            err = (abs(g_l_info)+1)*(1e10);
            return;
        end

    case 'hat'
        %         div = numel(xgrid_to_gl_vector)*2/numel(g_l_info);
        %
        %         for n = 1:numel(g_l_info)
        %             indices = xgrid_to_gl_vector==n;
        %             peak = max(x_grid(indices));
        %             half_width = max(x_grid(indices))-min(x_grid(indices));
        %             hat_vals = 1 - abs(peak - x_grid)/half_width;
        %             g_l_test = g_l_test + g_l_info(n)*hat_vals(:);
        %         end



        nhats = numel(g_l_info);
        intvl = round(numel(xgrid_to_gl_vector)/(nhats-1))*2;
        g_l_test = zeros((nhats+1)*intvl, 1);
        peak = (intvl+1)/2;
        half_width = (intvl-1)/2;
        hat_vals = 1 - abs(peak - (1:intvl)')/half_width;

        %first hat (half; falls off left boundary)
        g_l_test(1:intvl/2) = ...
            g_l_test(1:intvl/2) + g_l_info(1)*hat_vals(intvl/2+1:intvl);

        %add on all other hats
        for n = 0:nhats-2
            g_l_test(intvl*n/2+1:intvl/2*(n+2)) = ...
                g_l_test(intvl*n/2+1:intvl/2*(n+2)) + g_l_info(n+2)*hat_vals;
        end
        g_l_test = g_l_test(1:numel(x_grid));

        if sum(g_l_test<0)
            err = 1e15;
            return;
        end

    case 'polyn'
        dx = x_grid(2)-x_grid(1);
        l = x_grid(end);
        g_l_test = polyval(g_l_info, x_grid/l);

        if sum(g_l_test<=0)
            err = 1e10;
            return;
        end
        
    case 'sigmoid'
        l = x_grid(end);
        a = g_l_info(1); b = g_l_info(2); h = g_l_info(3); s = g_l_info(4);
        g_l_test = b + (a-b)./(1+exp((h-x_grid*100/l)/s));
        
        if sum(g_l_test<=0)
            err = 1e10;
            return;
        end
end

res_v = fhandle(g_l_test);

err = sum(sum((res_v(1,:,:)-v_tru(1,:,:)).^2,2),3);

subplot 131
%plot(reshape([res_v(1,:,:) v_tru(1,:,:)],size(v_tru,2),size(v_tru,3)*2,1))
plot([v_tru(1,:,1); res_v(1,:,1)]');
legend('true v', 'proposed v')
title('Comparison at x=0,\theta=\theta_1')
a = axis;
axis([a(1:2), -10 10])

subplot 132
semilogy(reshape((res_v(1,:,:)-v_tru(1,:,:)).^2,size(v_tru,2),size(v_tru,3)));
title(['Square error: ' num2str(err)], 'FontSize', 15)
a = axis;
axis([a(1:2), 1e-9 1e2])

subplot 133
plot([g_l_tru(:) g_l_test(:)]);
a = axis;
axis([a(1:2) 0 2*max(g_l_tru)])
switch basis_fcn
    case 'polyn'
        title({'Test g_L';char(poly2sym(round(g_l_info*4)/4))})
    case 'sigmoid'
        title({'Test g_L'; ['[a b h s]= ' num2str(g_l_info)]})
end
legend('True', 'Proposed')
drawnow;
end



function v = cab_eq(x_grid, dx, g_l, Nt, dt, G, stim, a, tau)
%function v = cab_eq(x_grid, g_l, G)
% solves the cable equation using a trapezoid scheme
%
%input: x_grid - points of dicretization in space
%       g_l - conductance at points in x_grid
%       Nt - number of time steps
%       G - temperature dependent cytoplasmic conductance (opposite of
%       resistivity)
%       stim - stimulus vector corresponding to time, given to the first
%           space component
%       a - cable radius
%       tau - time constant


Nx = numel(x_grid);
Ntemp = numel(G);
v = zeros(Nx,Nt,Ntemp); %voltage, mV

%solve forward for v using g5
e = ones(Nx,1);
S = spdiags([e -2*e e], -1:1, Nx, Nx)/dx/dx;
S(1,1) = -1/dx/dx;
S(Nx,Nx) = -1/dx/dx;

persistent n_eval
if isempty(n_eval)
    n_eval = 0;
else
    n_eval = n_eval + 1;
end

f_mod_no = 100;
if mod(n_eval,f_mod_no)==0
  fprintf('Function evaluation %d\n', n_eval);
end
for T = 1:Ntemp
    %fprintf('Current temp. No: %d of %d\n', T, Ntemp);

    %backward euler
    lambda_sq = a*G(T)/2./g_l(:);
    replambda_sq = repmat(lambda_sq(:),1,Nx);
    
    reptau = repmat(tau(:),1,Nx);
    B =(S.*replambda_sq - speye(Nx))./reptau;  %-dt*S*Gvec(T)/C + diag(1 + dt*g_l/C);
    
    v(:,1, T) = 0; %E_L(T); %initial conditions

    [L,U] = lu(speye(Nx)-dt*B/2);

    r = v(:,1, T);
    r(1) = (stim(2) + stim(1));%*dt/2;

    for n = 2:Nt
        v(:,n,T) = U \ ( L \ r);%U \ (L \ (v(:,n-1,T) + force_piece));
        r = 2*v(:,n,T) - r;
        r(1) = r(1) + (stim(n) + stim(n+1));%;*dt/2;
    end

    %if you want to see the actual simuation:
    %    subplot(2,3,T)
    %    plot3(tg,xg,v(:,plot_t_trunc_vec,T), 'b')
    %    xlabel('time')
    %    ylabel('space')
    %axis([0 t_final 0 l min(E_L)-5  max(E_L)+5]);

    %calculate V
    %V(:,T) = sum(v(:,:,T),2)*dt;
end

end
