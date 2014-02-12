%dtree1.m
%Ryan George
%VIGRE Summer 2010
%
%dtree1 simulates a dendritic tree of two branches as current is
%input in into one of the two branches.
%

loc = [1, 251, 750]; %values from scdend; correspond to first
                     %branch end, second branch end, soma

stim = struct('t1',1, 't2',2, 'amp',1e-4,'Tfin',10);
dt = .001;

if 1
for n = 1:3
    stim.loc = loc(n);
    [v_soma(:,n) tq teigs mvec] = scdend(dt, stim, 10);
end
end

ncomp = [250 250 250]; %# compartments per branch
Nt = ceil(stim.Tfin/dt);
tgrid = (1:Nt)*dt;
Neig = 12;

%find true eigenvectors
l = .025;
a = 1e-4;
h = .125;
R = .3;
g_L = 1/15;
lamsq = a/2/R/g_L;
L = l/sqrt(lamsq);
old_teigs = zeros(Neig,1);
old_teigs(1:2:Neig) = -((1:2:Neig).^2)*(pi^2)/4/L/L;
curguess = 2;
for n =2:2:Neig
    while ~old_teigs(n)
        curfind = fzero(@(z) (1+(z/.125)*tan(z*L))*tan(z*L)*(2*(a^(3/2)))+ ...
                        (a^(3/2))*(tan(z*L)-z/.125),curguess);
        if n==2
            old_teigs(n) = -curfind^2;
        elseif (curfind^2 + old_teigs(n-2))>1e-2
            old_teigs(n) = -curfind^2;
        end
        curguess = curguess+.5;
    end
end

old_tq = zeros(Neig, 3);
old_tq(1:2:Neig,1) = 1/sqrt(L); %a^(-3/4)/sqrt(L);
old_tq(1:2:Neig,2) = -1/sqrt(L); %a^(-3/4)/sqrt(L);
denom = (1+(sqrt(-old_teigs(2:2:Neig))/h).*tan(sqrt(-old_teigs(2:2:Neig))*L));
b_odd =((L +(L+1/h)./denom).^(-1/2)); %a^(-3/4)*((L +(L+1/h)./denom).^(-1/2));
old_tq(2:2:Neig,1:2) = repmat(b_odd,1,2);
old_tq(2:2:Neig, 3) = b_odd./denom;


th_0 = -(1:Neig)'; %((1:Neig)'.^2)*(pi^2)/4/(.25^2);
q_0 = ones(Neig,1);
param_0 = [q_0; th_0];
fits = somastim_modelfit(param_0, tgrid);
clear x fits
[x(:,3),resnorm(3)] = lsqcurvefit(@somastim_modelfit ,param_0,tgrid,v_soma(:,3));
fits(:,3) = somastim_modelfit(x(:,3), tgrid);
subplot(1,3,3)
plot(tgrid, [fits(:,3) v_soma(:,3)])

param0(Neig+1:2*Neig) = x(3,Neig+1:2*Neig);

for n = 1:2
    [x(:,n),resnorm(n)] = lsqcurvefit(@(p,xd) somastim_modelfit(p,xd,x(1:Neig,3)),...
                                      param_0,tgrid,v_soma(:,n));
    fits(:,n) = somastim_modelfit(x(:,n), tgrid);
    subplot(1,3,n)
    plot(tgrid, [fits(:,n) v_soma(:,n)])
end

%ezplot('(1+(z/.125)*tan(z*.25))*tan(.25*z)*2+tan(.25*z)-z/.125')