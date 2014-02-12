function comp_g_l = fourth_tempx()
% Ryan George
% VIGRE, Rice University
% This code naively recovers time-integrated V from initial and Neumann
% data, then seeks to replicate g_l(x) using this data
%

a = 1; %micro_m
a = a/(1e4); %cm
C_m = 1; %microF/cm^2

l = .1; %cm
dx = 5e-4; %cm
x_grid = (0:dx:l)';
Nx = numel(x_grid);

t_final = 50; %ms
dt = 5e-3; %cm
time = 0:dt:t_final;
Nt = numel(time);
mintemp = 32; %Celsius
dtemp = .25;
Ntemp = 20;
Tgrid = mintemp:dtemp:(mintemp+dtemp*(Ntemp-1));

g_0 = .3; %mS/cm^2
g_l = x_grid*0 + g_0; %must be a vector of size(x_grid)
proc_g_l = g_l*2*pi*a*dx;

%for Na:
%E_L = -40 + (0:Ntemp)';

%assume G linear
aff_params = fit_G(); %KOhms
G = @(T) aff_params(1)*T +aff_params(2); %KOhms
Gprime = @(T) aff_params(1);

Gvec = G(Tgrid)*pi*(a^2)/dx;

C = 2*pi*a*dx*C_m;
tau = C_m./g_l;

v = zeros(Nx,Nt,Ntemp); %voltage, mV

V = zeros(Nx,Ntemp); %integrated over time

stim = zeros(Nt+1,1); %current fed into soma end of cell
t1 = 1;
t2 = 2;
stim_indices = round(t1/dt): round(t2/dt);
stim_amp = 1e-4; %microV
mfac = (dt/2)/(2*pi*a*dx)/C_m;%(dt/2)/(2*pi*a*dx)/C_m;
stim(stim_indices) = stim_amp*mfac;
%stim([round(t1/dt)-1,round(t2/dt)+1])=  stim_amp*mfac/2;
%stim = stim/(1e3); %convert to mV;

%solve forward for v using g5
e = ones(Nx,1);
S = spdiags([e -2*e e], -1:1, Nx, Nx)/dx/dx;
S(1,1) = -1/dx/dx;
S(Nx,Nx) = -1/dx/dx;

clf;
title('Forward-Constructed V', 'FontSize', 15)
mod_num = 30;
plot_t_trunc_vec = 1:mod_num:Nt;
[tg xg] = meshgrid(time(plot_t_trunc_vec), x_grid);

%e1 = zeros(Nx,1); e1(1)=1;



for T = 1:Ntemp
    T
    %backward euler
    lambda_sq = a*G(Tgrid(T))/2./g_l;
    replambda_sq = repmat(lambda_sq(:),1,Nx);
    B =(S.*replambda_sq - eye(Nx))/tau(T);  %-dt*S*Gvec(T)/C + diag(1 + dt*g_l/C);
    v(:,1, T) = 0; %E_L(T); %initial conditions

    [L,U] = lu(eye(Nx)-dt*B/2);

    r = v(:,1, T);
    r(1) = (stim(2) + stim(1));%*dt/2;

    for n = 2:Nt
        v(:,n,T) = U \ ( L \ r);%U \ (L \ (v(:,n-1,T) + force_piece));
        r = 2*v(:,n,T) - r;
        r(1) = r(1) + (stim(n) + stim(n+1));%;*dt/2;
        if n*dt>1
        end
    end

    %if you want to see the actual simuation:
%    subplot(2,3,T)
%    plot3(tg,xg,v(:,plot_t_trunc_vec,T), 'b')
%    xlabel('time')
%    ylabel('space')
    %axis([0 t_final 0 l min(E_L)-5  max(E_L)+5]);

    %calculate V
    V(:,T) = sum(v(:,:,T),2)*dt;
end

clear B L U replambda_sq


%subplot 121
%[xg Tg] = meshgrid(Tgrid, x_grid);
%surf(xg, Tg, V)
%a = axis;
%axis([mintemp max(Tgrid) 0 l 0 3])

%analytical solution
if 0
nterms = 40;
q = zeros(Nx,nterms);
q(:,1) = 1/sqrt(l);
q(:,2:end) = sqrt(2/l)*cos(x_grid*(1:(nterms-1))*pi/l);
q0_form = reshape(q(1,:),1,1,nterms);
q0rep = repmat(q0_form,Nx,Nt);
qrep = repmat(reshape(q, [Nx,1,nterms]),1,Nt);
q_terms = q0rep.*qrep;

clear q q0rep qrep

sc_theta = -((0:nterms-1)*pi/l).^2;

true_v = zeros(size(v));
Ctau = C_m/g_0;%for constant g_l


for T = 1:Ntemp
    Clambda_sq = a*G(Tgrid(T))/2/g_0; %for constant g_l
    sc_zeta = reshape((Clambda_sq*sc_theta - 1)/Ctau,1,1,nterms);

    stim_int = zeros(Nx,Nt,nterms); %stimulus integral
    %    stim_int(:,time<=t1) = 0;
    stim_on = find((time>t1).*(time<=t2));

    ib_time = time(stim_on);
    rep_ib_time = repmat(ib_time, [Nx 1 nterms]);
    repzeta = repmat(sc_zeta,Nx,numel(ib_time));
    stim_int(:,stim_on,:) = stim_amp*(exp(repzeta.*(rep_ib_time-t1))-1)./repzeta;

    post_stim = find(time>t2);
    ps_time = time(post_stim);
    rep_ps_time = repmat(ps_time, [Nx 1 nterms]);
    repzeta = repmat(sc_zeta,Nx,numel(ps_time));
    stim_int(:,post_stim,:) = stim_amp*(exp(repzeta.*(rep_ps_time-t1))-exp(repzeta.*(rep_ps_time-t2)))./repzeta;

    an_v(:,:,T) = sum(q_terms.*stim_int,3)/C_m/2/a/pi;
    plot3(tg,xg,an_v(:,plot_t_trunc_vec,T), 'b')
end

an_V = reshape(sum(an_v,2)*dt,[Nx,Ntemp]); %actually

end

%recontruct V from data
sV = zeros(Nx,Ntemp);
sV(1,:) = V(1,:);

I = stim_amp*(t2-t1);%sum(stim)*dt;

sV(2,:) = sV(1,:) - I./Gvec;

rhs = zeros(Ntemp,1);
cheat_rhs = zeros(Ntemp,1);
itv = Tgrid(2 : end-1); %inner temp vec





for n = 2:Nx-1

    cheat_mm1_diag = [-V(n,2:end-1).*G(itv)/2/dx/dx/dtemp   -V(n,Ntemp)*G(Tgrid(Ntemp))/dx/dx/dtemp];

    cheat_me1_diag = [V(n,1)*Gprime(mintemp)/dx/dx - G(mintemp)/dtemp/dx/dx*(V(n,1) + V(n,2) - V(n,1) )...
        V(n,2:end-1).*Gprime(itv)/dx/dx - G(itv)/2/dx/dx/dtemp.*(V(n,3:end)-V(n,1:end-2))...
        V(n,end)*Gprime(Tgrid(Ntemp))/dx/dx + G(Tgrid(Ntemp))/dtemp/dx/dx*(V(n,end) - V(n,end) + V(n,end-1))];

    cheat_mp1_diag = [V(n,1)*G(mintemp)/dtemp/dx/dx   V(n,2:end-1).*G(itv)/2/dtemp/dx/dx  ];

    cheat_A = diag(cheat_mm1_diag,-1) ... %for V(n+1,m-1)
        + diag(cheat_me1_diag)... %for V(n+1,m)
        +diag(cheat_mp1_diag,1); %for V(n+1,m+1)

    
    cheat_rhs(1) = G(Tgrid(1))/dx/dx/dtemp.*(-2*V(n,1)+V(n-1,1)).*(V(n,2)-V(n,1))...
        -V(n,1).*(G(Tgrid(1))/dx/dx/dtemp.*(-2*V(n,2)+V(n-1,2)+2*V(n,1)-V(n-1,1))...
                + Gprime(Tgrid(1))/dx/dx.*(-2*V(n,1)+V(n-1,1)));
    
    cheat_rhs(2:Ntemp-1) = G(itv)/dx/dx/2/dtemp.*(-2*V(n,2:end-1)+V(n-1,2:end-1)).*(V(n,3:end)-V(n,1:end-2))...
        -V(n,2:end-1).*(G(itv)/dx/dx/2/dtemp.*(-2*V(n,3:end)+V(n-1,3:end)+2*V(n,1:end-2)-V(n-1,1:end-2))...
                      + Gprime(itv)/dx/dx.*(-2*V(n,2:end-1)+V(n-1,2:end-1)));
            
    cheat_rhs(Ntemp) = G(Tgrid(end))/dx/dx/dtemp.*(-2*V(n,end)+V(n-1,end)).*(V(n,end)-V(n,end-1))...
        -V(n,end).*(G(Tgrid(end))/dx/dx/dtemp.*(-2*V(n,end)+V(n-1,end)+2*V(n,end-1)-V(n-1,end-1))...
                  + Gprime(Tgrid(end))/dx/dx.*(-2*V(n,end)+V(n-1,end)));
    
    
              
    mm1_diag = [-sV(n,2:end-1).*G(itv)/2/dx/dx/dtemp   -sV(n,Ntemp)*G(Tgrid(Ntemp))/dx/dx/dtemp];

    me1_diag = [sV(n,1)*Gprime(mintemp)/dx/dx - G(mintemp)/dtemp/dx/dx*(sV(n,1) + sV(n,2) - sV(n,1) )...
        sV(n,2:end-1).*Gprime(itv)/dx/dx - G(itv)/2/dx/dx/dtemp.*(sV(n,3:end)-sV(n,1:end-2))...
        sV(n,end)*Gprime(Tgrid(Ntemp))/dx/dx + G(Tgrid(Ntemp))/dtemp/dx/dx*(sV(n,end) - sV(n,end) + sV(n,end-1))];

    mp1_diag = [sV(n,1)*G(mintemp)/dtemp/dx/dx   sV(n,2:end-1).*G(itv)/2/dtemp/dx/dx  ];

    A = diag(mm1_diag,-1) ... %for sV(n+1,m-1)
        + diag(me1_diag)... %for sV(n+1,m)
        +diag(mp1_diag,1); %for sV(n+1,m+1)

    
    rhs(1) = G(Tgrid(1))/dx/dx/dtemp.*(-2*sV(n,1)+sV(n-1,1)).*(sV(n,2)-sV(n,1))...
        -sV(n,1).*(G(Tgrid(1))/dx/dx/dtemp.*(-2*sV(n,2)+sV(n-1,2)+2*sV(n,1)-sV(n-1,1))...
                + Gprime(Tgrid(1))/dx/dx.*(-2*sV(n,1)+sV(n-1,1)));
    
    rhs(2:Ntemp-1) = G(itv)/dx/dx/2/dtemp.*(-2*sV(n,2:end-1)+sV(n-1,2:end-1)).*(sV(n,3:end)-sV(n,1:end-2))...
        -sV(n,2:end-1).*(G(itv)/dx/dx/2/dtemp.*(-2*sV(n,3:end)+sV(n-1,3:end)+2*sV(n,1:end-2)-sV(n-1,1:end-2))...
                      + Gprime(itv)/dx/dx.*(-2*sV(n,2:end-1)+sV(n-1,2:end-1)));
            
    rhs(Ntemp) = G(Tgrid(end))/dx/dx/dtemp.*(-2*sV(n,end)+sV(n-1,end)).*(sV(n,end)-sV(n,end-1))...
        -sV(n,end).*(G(Tgrid(end))/dx/dx/dtemp.*(-2*sV(n,end)+sV(n-1,end)+2*sV(n,end-1)-sV(n-1,end-1))...
                  + Gprime(Tgrid(end))/dx/dx.*(-2*sV(n,end)+sV(n-1,end)));
                  
    matrix = A; 
%     if n < 100
%         divvec = cheat_rhs;
%     else    
         divvec = rhs;
%     end
    
    sV(n+1,:) = matrix \ divvec;


end
%subplot 122
%surf(xg, Tg, sV)

figure(3)
zlim = max(max(V))*2;
%look at recontruction
[VTg Vxg ] = meshgrid(Tgrid, x_grid);
subplot 121
plot3( Vxg, VTg,V)
title('Simluated V', 'FontSize', 15)
axis([0 l mintemp Tgrid(end) 0 zlim])
subplot 122
plot3( Vxg, VTg,sV)
title('Reconstructed V', 'FontSize', 15)
axis([0 l mintemp Tgrid(end) 0 zlim])

%reconstruct g - note, for now we use true V
Gmat = repmat(G(Tgrid(:))',Nx,1);
naive_rec_g_l =   mean((S*sV).*Gmat./sV,2)*a/2;
smooth_rec_g_l = ( naive_rec_g_l(1:end-1) + naive_rec_g_l(2:end))/2;
smooth_rec_g_l(Nx) = smooth_rec_g_l(Nx-1);


ornate_S =  spdiags([-e 16*e -30*e 16*e -e], -2:2, Nx, Nx)/dx/dx/12;
ornate_S(1:2,1:4) = [-1 1 0 0     %clean up the edges
                     1 -2 1 0]/dx/dx;
ornate_S(Nx-1:Nx,Nx-3:Nx) = [0 1 -2 1
                             0 0 1 -1 ]/dx/dx;

second_estimate_rec_g_l = mean((ornate_S*sV).*Gmat./sV,2)*a/2;

smooth_se_rec_g_l = (second_estimate_rec_g_l(1:end-1) + second_estimate_rec_g_l(2:end))/2;
smooth_se_rec_g_l(Nx) = smooth_se_rec_g_l(Nx-1);

comp_g_l = [g_l naive_rec_g_l smooth_rec_g_l second_estimate_rec_g_l smooth_se_rec_g_l];

figure(34)
plot(x_grid,comp_g_l)
legend('true', 'recovered', 'smoothed', 'better \partial_{x}^{2} estimate', 'better estimate smoothed')
axis([0 l 0 1])

figure(35)
plot(x_grid,comp_g_l(:,[1 3]))
legend('true', 'best recovered estimate')
description = ['g_l recovered: dx = ' num2str(dx) ', dt = ' num2str(dt), ', Ntemp = ' num2str(Ntemp) ];
title(description, 'FontSize', 14)
axis([0 l 0 max(g_l)*2])
end