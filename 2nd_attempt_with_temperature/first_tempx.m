function first_tempx
% Ryan George
% VIGRE, Rice University
% This code naively recovers time-integrated V from initial and Neumann 
% data, then seeks to replicate g_l(x) using this data
%

a = 1; %micro_m
a = a/(1e4); %cm
C_m = 1; %microF/cm^2

l = .1; %cm
dx = 1e-2; %cm
x_grid = (0:dx:l)';
Nx = numel(x_grid);

t_final = 2; %ms
dt = 5e-3; %cm
time = 0:dt:t_final;
Nt = numel(time);
mintemp = 32; %Celsius
dtemp = 3;
Ntemp = 6;
Tgrid = mintemp:dtemp:(mintemp+dtemp*(Ntemp-1));

g_0 = 1/15; %mS/cm^2
pre_g_l = x_grid*0 + g_0; %must be a vector of size(x_grid)
g_l = pre_g_l*2*pi*a*dx;

%for Na:
E_L = -40 + (0:Ntemp)';

%assume G linear
aff_params = fit_G(); %KOhms
G = @(temp) aff_params(1)*temp +aff_params(2); %KOhms
Gprime = @(temp) aff_params(1);

C = 2*pi*a*dx*C_m;
Gvec = G(Tgrid)*(pi*a^2)^2/dx;


v = zeros(Nx,Nt,Ntemp); %voltage, mV

V = zeros(Nx,Ntemp); %integrated over time

stim = zeros(Nt,1); %current fed into soma end of cell
stim(20:40) = 60; %mV

%solve forward for v using g5
e = ones(Nx,1);
S = spdiags([e -2*e e], -1:1, Nx, Nx)/dx/dx;
S(1,1) = -1/dx/dx;
S(Nx,Nx) = -1/dx/dx;

figure(1); clf;
title('Forward-Constructed V', 'FontSize', 15)

for T = 1:Ntemp
    %backward euler
    D =  -dt*S*Gvec(T) + diag(C + dt*g_l);
    v(:,1, T) = E_L(T); %initial conditions
    
    force_piece = dt*g_l*E_L(T);

    for n = 2:Nt
        force_piece(1) = dt*g_l(1)*E_L(T) + dt*stim(n); %stim should be scaled
        v(:,n,T) = D \ (C*v(:,n-1,T) + force_piece);
    end

    %if you want to see the actual simuation:
    mod_num = 3;
    plot_t_trunc_vec = 1:mod_num:Nt;
    [tg xg] = meshgrid(time(plot_t_trunc_vec), x_grid);
    subplot(2,3,T)
    plot3(tg,xg,v(:,plot_t_trunc_vec,T))
    xlabel('time')
    ylabel('space')
    axis([0 t_final 0 l min(E_L)-5  max(E_L)+5]);

    %calculate V
    V(:,T) = sum(v(:,:,T)-E_L(T),2)*dt;
end

%subplot 121
%[xg Tg] = meshgrid(Tgrid, x_grid);
%surf(xg, Tg, V)
%a = axis;
%axis([mintemp max(Tgrid) 0 l 0 3])


%recontruct V from data
sV = zeros(Nx,Ntemp);
sV(1,:) = V(1,:);

I = sum(stim)*dt;

sV(2,:) = sV(1,:) + dx * I;

rhs = zeros(Ntemp,1);
itv = Tgrid(2 : end-1); %inner temp vec

for n = 2:Nx-1
    
    mm1_diag = [sV(n,2:end-1).*G(itv)/2/dx/dx/dtemp   -sV(n,Ntemp)*G(Tgrid(Ntemp))/dx/dx/dtemp];
    me1_diag = [sV(n,1)*Gprime(mintemp)/dx/dx - G(mintemp)/dtemp/dx/dx*(sV(n,1) + sV(n,2) - sV(n,1) )...
                sV(n,2:end-1).*Gprime(itv)/dx/dx - G(itv)/2/dx/dx/dtemp.*(sV(n,3:end)-sV(n,1:end-2))...
                sV(n,end)*Gprime(Tgrid(Ntemp))/dx/dx + G(mintemp)/dtemp/dx/dx*(sV(n,end) - sV(n,end) + sV(n,end-1))];
    mp1_diag = [sV(n,1)*G(mintemp)/dtemp/dx/dx   sV(n,2:end-1).*G(itv)/2/dtemp/dx/dx  ];
    
    A = diag(mm1_diag,-1) ... %for sV(n+1,m-1)
        + diag(me1_diag)... %for sV(n+1,m)
        +diag(mp1_diag,1); %for sV(n+1,m+1)
    
    
    rhs(1) = sV(n,1)*(Gprime(mintemp)/dx/dx*(2*sV(n,1) + sV(n-1,1))...
        + G(mintemp)/dtemp/dx/dx*(2*sV(n,2)-sV(n-1,2)-2*sV(n,1)+sV(n-1,1) ) )...
        + G(mintemp)/dtemp/dx/dx*(-2*sV(n,1)+sV(n-1,1))*(sV(n,2)-sV(n,1));
        
    rhs(2:Ntemp-1) = (sV(n,2:end-1)/dx/dx.*(Gprime(itv).*(2*sV(n,2:end-1)-sV(n-1,2:end-1))...
        + G(itv)/dtemp.*(2*sV(n,3:end)-sV(n-1,3:end)-2*sV(n,1:end-2)+sV(n-1,1:end-2))) ...
        +G(itv)/dx/dx/dtemp.*(-2*sV(n,2:end-1)+sV(n-1,2:end-1)).*(sV(n,3:end)-sV(n,1:end-2)))';
    
    rhs(Ntemp) = sV(n,end)*(Gprime(Tgrid(Ntemp))/dx/dx*(2*sV(n,end) + sV(n-1,end))...
        + G(mintemp)/dtemp/dx/dx*(2*sV(n,end)-sV(n-1,end)-2*sV(n,end-1)+sV(n-1,end-1) ) )...
        + G(mintemp)/dtemp/dx/dx*(-2*sV(n,end)+sV(n-1,end))*(sV(n,end)-sV(n,end-1));
    
    sV(n+1,:) = A \ rhs;
        
    
end
%subplot 122
%surf(xg, Tg, sV)

end