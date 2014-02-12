function rec_g_l = fourth_tempx()
% Ryan George
% VIGRE, Rice University
% This code naively recovers time-integrated V from initial and Neumann 
% data, then seeks to replicate g_l(x) using this data
%

a = 1; %micro_m
a = a/(1e4); %cm
C_m = 1; %microF/cm^2
%C_m = C_m/(1e4) %mF/cm^2

l = .1; %cm
dx = 1e-3; %cm
x_grid = (0:dx:l)';
Nx = numel(x_grid);

t_final = 10; %ms
dt = 5e-3; %cm
time = 0:dt:t_final;
Nt = numel(time);
mintemp = 32; %Celsius
dtemp = 3;
Ntemp = 6;
Tgrid = mintemp:dtemp:(mintemp+dtemp*(Ntemp-1));

g_0 = 1/15; %mS/cm^2
g_l = x_grid*10 + g_0; %must be a vector of size(x_grid)
proc_g_l = g_l*2*pi*a*dx;

%for Na:
%E_L = -40 + (0:Ntemp)';

%assume G linear
aff_params = fit_G(); %KOhms
G = @(T) aff_params(1)*T +aff_params(2); %KOhms
Gprime = @(T) aff_params(1);

Gvec = G(Tgrid)*pi*a^2/dx;

C = 2*pi*a*dx*C_m;
tau = C_m./g_l;

v = zeros(Nx,Nt,Ntemp); %voltage, mV

V = zeros(Nx,Ntemp); %integrated over time

stim = zeros(Nt,1); %current fed into soma end of cell
stim_on = 1;
stim_off = 2;
stim_indices = round(stim_on/dt): round(stim_off/dt);
stim_amp = 1e-2; %microV
stim(stim_indices) = stim_amp; 
stim([round(stim_on/dt)-1,round(stim_off/dt)+1])=  stim_amp/2;
stim = stim/(1e3); %convert to mV;

%solve forward for v using g5
e = ones(Nx,1);
S = spdiags([e -2*e e], -1:1, Nx, Nx)/dx/dx;
S(1,1) = -1/dx/dx;
S(Nx,Nx) = -1/dx/dx;

figure(1); clf;
title('Forward-Constructed V', 'FontSize', 15)  
e1 = zeros(Nx,1); e1(1)=1;

for T = 1:Ntemp
    %backward euler
    lambda_sq = a*G(T)/2./g_l;
    replambda_sq = repmat(lambda_sq(:),1,Nx);
    B =(S.*replambda_sq - eye(Nx))/tau(T);  %-dt*S*Gvec(T)/C + diag(1 + dt*g_l/C);
    v(:,1, T) = 0; %E_L(T); %initial conditions
    
    [L,U] = lu(eye(Nx)-dt*B);
    
    for n = 2:Nt
        v(:,n,T) = U \ ( L \ (v(:,n-1,T) + e1*dt* stim(n)/C));%U \ (L \ (v(:,n-1,T) + force_piece));
    end

    %if you want to see the actual simuation:
    mod_num = 30;
    plot_t_trunc_vec = 1:mod_num:Nt;
    [tg xg] = meshgrid(time(plot_t_trunc_vec), x_grid);
    subplot(2,3,T)
    plot3(tg,xg,v(:,plot_t_trunc_vec,T), 'b')
    xlabel('time')
    ylabel('space')
    %axis([0 t_final 0 l min(E_L)-5  max(E_L)+5]);

    %calculate V
    V(:,T) = sum(v(:,:,T),2)*dt;
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

%reconstruct g - note, for now we use true V
Gmat = repmat(Gvec(:)',Nx,1);
rec_g_l = [g_l  mean((S*V).*Gmat./V,2)];

end