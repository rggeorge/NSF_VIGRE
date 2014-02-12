%somastim_modelfit.m
%calculates the expected fit from equation 2 of Dendrite_g_L.pdf

function out = somastim_modelfit(params, tgrid, q_obs)

t1 = 1;
t2 = 2;
g_L = 1/15;
tau = 1/g_L;
A_s = 4*pi*1e-6;
a = 1e-4;
R2 = .3;
lam = a/2/R2/g_L;
h = a*lam*2*pi/A_s;
subh = a*a*pi/R2;

Npar = numel(params)/2;
q =  params(1:Npar); q = q(:);
th = params((Npar+1):end); th = th(:);

Nt = numel(tgrid);
dt = tgrid(2)-tgrid(1);
Tgrid = tgrid/tau;

%calculate integral (of stimulus)
intgl = zeros(Npar, Nt);
onind =(t1/dt):(t2/dt);
intgl(:,onind) = 1 - exp(th*(Tgrid(onind)-t1/tau));

it = th*(Tgrid(onind(end)+1:end)-t2/tau);
it2 = th*(Tgrid(onind(end)+1:end)-t1/tau);
intgl(:,t2/dt+1:Nt) = exp(it) - exp(it2);

%intgl(:,t2/dt+1:Nt) = exp(th*(Tgrid(onind(end)+1:end)-t2/tau))...
% - exp(th*(Tgrid(onind(end)+1:end)-t1/tau));

amp = 1e-4;%amps*10e-8

if nargin>2
    out = intgl'*(q.*q_obs./(1-th))*(1/subh)*(1e-4);
else
    out = intgl'*((q.^2)./(1-th))*(1/subh)*amp;
end