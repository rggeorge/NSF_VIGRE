function  [g_L_output error] = bf2(basis_fcn, Ntemp, Trange, Nbf, bf_mod_no)
% Ryan George
% VIGRE, Rice University
% This code naively recovers time-integrated V from initial and Neumann
% data, then seeks to replicate g_L(x) using this data
%
% duplicted from bf1.m on 6/15
% goals: emulate Cox (Adjoint Method)
%        find effect of adding more temperatures

if nargin < 1
    basis_fcn = 'step';
end

if nargin<2
   Ntemp = 1;
end

if nargin<3
    Trange = 1;
end

if nargin < 4
    Nbf = 10;
end

if nargin < 5
    bf_mod_no = 4; %number of steps in g_L step function
end


a = 1; %micro_m
a = a/(1e3); %cm
C_m = 1; %microF/cm^2
l = .1; %cm
        %dx = 2.5e-3; %cm

Nx = Nbf*bf_mod_no;
x_grid = linspace(0,l,Nx)';
dx = l/Nx;

t_final = 6; %ms
dt = 2e-3; %cm from Fohlmeister: times digitized at 10, 20, 40 kHz 
  %                                            dt = .1, .05, .025 ms
time = 0:dt:t_final;
Nt = numel(time);

maxtemp = 32;
mintemp = 27;

if nargin>1
    mintemp = maxtemp - Trange;
end

Tgrid = linspace(mintemp,maxtemp, Ntemp);
%dtemp = Tgrid(2)-Tgrid(1);

g_L_fcn = @(an_x_grid) ((an_x_grid).^2)*1e3 + .5;
g_L = g_L_fcn(x_grid); %must be a vector of size(x_grid)

C = 2*pi*a*dx*C_m;
tau = C_m./g_L;

stim_amp = 1/800;%1e-4; %microV
mfac = (dt/2)/(2*pi*a*dx)/C_m;%(dt/2)/(2*pi*a*dx)/C_m;
%stim = .3*max(time-1,0).*exp(-max(time-1,0)/2);
stim = (time>1).*(time<5);

stim(end+1)=0; %to avoid an error in the last loop
stim = stim*stim_amp*mfac;

cab = struct("x_grid", x_grid, "dx", dx, "Nt", ...
             Nt, "dt", dt, "stim", stim, "C_m", C_m, "Tgrid", Tgrid);
runv = @(test_g_L) cab_eq(x_grid, dx, test_g_L, Nt, dt, ...
                          stim, a, C_m, Tgrid);

v_tru = runv(g_L);

if 1
clf;
title('Forward-Constructed V', 'FontSize', 15)
mod_num = 30;
plot_t_trunc_vec = 1:mod_num:Nt;
[tg xg] = meshgrid(time(plot_t_trunc_vec), x_grid);
plot3(tg,xg,v_tru(:,plot_t_trunc_vec,1), 'b')
end



%test and recover g_L

switch basis_fcn
    case 'step'
        Nbf = ceil(Nx/bf_mod_no);        
    case 'hat'
        Nbf = ceil(Nx/bf_mod_no) + 2;
    case 'polyn'
        Nbf = bf_mod_no;
end

minfcn = @(g_L_info) g_L_minimization(g_L_info, basis_fcn, bf_mod_no, x_grid, v_tru, runv, g_L);

options = optimset('MaxFunEvals', 1e3, 'TolFun', 1e-16, 'Display', 'final');

tic
figure(2)
switch basis_fcn

    case 'step'
        initial_guess = ones(Nbf,1)*mean(g_L).*(1+.05*randn(Nbf,1));
        g_L_rec = fminunc(@(g_L_info) minfcn(g_L_info), initial_guess, options);

    case 'hat'
        initial_guess = ones(Nbf,1)*g_0+.05;
        %       initial_guess = g_L_fcn(linspace(0,l,Nbf))*1.35;
 
%         g_L_rec = fmincon(@(g_L_info) minfcn(g_L_info), initial_guess,[],[],[],[],...
%             ones(Nbf,1)*.05,ones(Nbf,1)*3*max(g_L),[], options);
        g_L_rec = fminunc(@(g_L_info) minfcn(g_L_info), initial_guess, options);
    case 'polyn'
        initial_guess = [zeros(Nbf-1,1); 1];
        initial_guess = g_L(linspace(0,l,Nbf))+.1;

        g_L_rec = fminunc(@(g_L_info) minfcn(g_L_info), initial_guess);

    case 'sigmoid'
        initial_guess = [1 1 1 1]; %[a b h s]
        g_L_rec = fminunc(@(g_L_info) minfcn(g_L_info), initial_guess);
end

elapsed_time = toc

figure(4); clf;
%center_pts = (xgl_pts(2:end) + xgl_pts(1:end-1))/2;



switch basis_fcn
    case 'step'
        xstep = zeros(Nbf*2, 1);
        step_g_L = zeros(Nbf*2, 1);
        for n = 1:Nbf
            xstep(n*2:n*2+1) = n*bf_mod_no*dx;
            step_g_L(n*2-1:n*2) = g_L_rec(n);
            initial_g_L((n*2-1):n*2) = initial_guess(n);
        end

        plot(x_grid, g_L, xstep(1:Nbf*2), step_g_L, 'g', xstep(1:Nbf*2), initial_g_L, 'c-o')

        for n = 1:Nbf-1
            g_L_fin( ((n-1)*bf_mod_no+1):n*bf_mod_no ) = g_L_rec(n);
        end
            g_L_fin( ((Nbf-1)*bf_mod_no+1) : Nx ) = g_L_rec(Nbf);
            %^           g_L_test( ((Nbf-1)*bf_mod_no+1:numel(x_grid))=g_L_info(Nbf);



            
    case 'hat'
        intvl = round(Nx/(Nbf-1))*2;
        g_L_fin = zeros((Nbf+1)*intvl, 1);
        peak = (intvl+1)/2;
        half_width = (intvl-1)/2;
        hat_vals = 1 - abs(peak - (1:intvl)')/half_width;

        %first hat (half; falls off left boundary)
        g_L_fin(1:intvl/2) = ...
            g_L_fin(1:intvl/2) + g_L_rec(1)*hat_vals(intvl/2+1:intvl);
        

        %add on all other hats
        for n = 0:Nbf-2
            g_L_fin(intvl*n/2+1:intvl/2*(n+2)) = ...
                g_L_fin(intvl*n/2+1:intvl/2*(n+2)) + g_L_rec(n+2)*hat_vals;
        end
        g_L_fin = g_L_fin(1:Nx);
        plot(x_grid, [g_L(:)  g_L_fin(:)]);

    case 'polyn'
        g_L_fin = polyval(g_L_rec, x_grid/l);
        plot(x_grid, [g_L(:)  g_L_fin(:)]);

    case 'sigmoid'
        a = g_L_rec(1); b = g_L_rec(2); h = g_L_rec(3); s = g_L_rec(4);
        g_L_fin = b + (a-b)./(1+exp((h-x_grid*100/l)/s));
        plot(x_grid, [g_L(:)  g_L_fin(:)]);
        
end


title({'Recovering Conductance: Brute Force on Lateral Overdetermination';...
    ['dx = ' num2str(dx) ', dt = ' num2str(dt) ', #' basis_fcn ...
    's = ' num2str(Nbf)];...
    ['\theta range = [' num2str(Tgrid(1)) ' ' num2str(Tgrid(end)) ']'...
    ', #\theta = ' num2str(Ntemp)]}, 'FontSize', 15)
axis([0 l 0 2*max(g_L)])
xlabel('Length (cm)', 'FontSize', 12)
ylabel('g_L (mS/cm^2)', 'FontSize', 12)
legend('true g_L', 'recovered g_L');


error = sum((g_L_fin(:) - g_L(:)).^2);
g_L_output = g_L_rec(:);

end




function err = g_L_minimization(g_L_info, basis_fcn, bf_mod_no, x_grid, v_tru, fhandle, g_L_tru)


switch basis_fcn


    case 'step'
        Nbf = numel(g_L_info);        
        g_L_test = zeros(size(x_grid));
        for n = 1:Nbf-1
            g_L_test( ((n-1)*bf_mod_no+1):n*bf_mod_no  ) = g_L_info(n);
        end
        g_L_test((Nbf-1)*bf_mod_no+1:numel(x_grid))=g_L_info(Nbf);
        
        %        if sum(g_L_info<=0)
        %    err = (abs(g_L_info)+1)*(1e10);
        %    return;
        %end

    case 'hat'

        nhats = numel(g_L_info);
        intvl = round(numel(xgrid_to_gl_vector)/(nhats-1))*2;
        g_L_test = zeros((nhats+1)*intvl, 1);
        peak = (intvl+1)/2;
        half_width = (intvl-1)/2;
        hat_vals = 1 - abs(peak - (1:intvl)')/half_width;

        %first hat (half; falls off left boundary)
        g_L_test(1:intvl/2) = ...
            g_L_test(1:intvl/2) + g_L_info(1)*hat_vals(intvl/2+1:intvl);

        %add on all other hats
        for n = 0:nhats-2
            g_L_test(intvl*n/2+1:intvl/2*(n+2)) = ...
                g_L_test(intvl*n/2+1:intvl/2*(n+2)) + g_L_info(n+2)*hat_vals;
        end
        g_L_test = g_L_test(1:numel(x_grid));

        %        if sum(g_L_test<0)
        %    err = 1e15;
        %    return;
        %end

    case 'polyn'
        dx = x_grid(2)-x_grid(1);
        l = x_grid(end);
        g_L_test = polyval(g_L_info, x_grid/l);

        if sum(g_L_test<=0)
            err = 1e10;
            return;
        end
        
    case 'sigmoid'
        l = x_grid(end);
        a = g_L_info(1); b = g_L_info(2); h = g_L_info(3); s = g_L_info(4);
        g_L_test = b + (a-b)./(1+exp((h-x_grid*100/l)/s));
        
        if sum(g_L_test<=0)
            err = 1e10;
            return;
        end
end

res_v = fhandle(g_L_test);

x1 = .075;
[minval x1_index] = min(abs(x_grid-x1));

err = sum(sum((res_v(1,:,:)-v_tru(1,:,:)).^2,2),3);%+...
%      sum(sum((res_v(x1_index,:,:)-v_tru(x1_index,:,:)).^2,2),3);


f_mod_no = 20;
if mod(n_eval,f_mod_no)==0
    fprintf('Function evaluation %d\n', n_eval);
    subplot 131
    plot([v_tru(1,:,1); res_v(1,:,1)]');
    legend('true v', 'proposed v')
    title('Comparison at x=0,\theta=\theta_1')
    a = axis;

    subplot 132
    semilogy(reshape((res_v(1,:,:)-v_tru(1,:,:)).^2,size(v_tru,2),size(v_tru,3)));
    title(['Square error: ' num2str(err)], 'FontSize', 15)
    a = axis;
    axis([a(1:2), 1e-9 1e2])

    subplot 133
    plot([g_L_tru(:) g_L_test(:)]);
    a = axis;
    axis([a(1:2) 0 2*max(g_L_tru)])
    switch basis_fcn
      case 'polyn'
          title({'Test g_L';char(poly2sym(round(g_L_info*4)/4))})
      case 'sigmoid'
          title({'Test g_L'; ['[a b h s]= ' num2str(g_L_info)]})
    end
legend('True', 'Proposed')
drawnow;
end

end


function v = cab_eq(x_grid, dx, g_L, Nt, dt, stim, a, C_m, Tgrid)
%function v = cab_eq(x_grid, g_L, G)
% solves the cable equation using a trapezoid scheme
%
%input: x_grid - points of dicretization in space
%       g_L - conductance at points in x_grid
%       Nt - number of time steps
%       G - temperature dependent cytoplasmic conductance (opposite of
%       resistivity)
%       stim - stimulus vector corresponding to time, given to the first
%           space component
%       a - cable radius
%       tau - time constant


Nx = numel(x_grid);
Ntemp = numel(Tgrid);
v = zeros(Nx,Nt,Ntemp); %voltage, mV

[all_ts_g_L G] = apply_q10(g_L, Tgrid);  %[mS/cm^2   Kohm*cm]

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

Nx = numel(x_grid);

for T = 1:Ntemp
    %temperature-specific g_L
    ts_g_L = all_ts_g_L(:,T);

    B =(S.*a*G(T)/2 - spdiags(ts_g_L,0,Nx,Nx))/C_m;
    
    v(:,1, T) = 0; %E_L(T); %initial conditions

    [L,U] = lu(speye(Nx)-dt*B/2);
    
    r = v(:,1, T);
    r(1) = (stim(2) + stim(1));  %premultiplied by dt/2;
    
    for n = 2:Nt
        v(:,n,T) = U \ ( L \ r);%U \ (L \ (v(:,n-1,T) + force_piece));
        r = 2*v(:,n,T) - r;
        r(1) = r(1) + (stim(n) + stim(n+1));  %premultiplied by dt/2;
    end

end

end
