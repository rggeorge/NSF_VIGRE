% Ryan George
% CAAM VIGRE, Rice University
% Summer 2010
%
% qsolve.m solves for the eigenfunctions given by the differential
% equation q''(x) + g_L(x)q(x) + th * q(x) = 0

function [tq teigs] = qsolve(numreq)

if nargin < 1
    numreq = 8;
end

a = 1e-4*[1 1 1];
R_a = 0.3;
ell = [2.5 2.5 2.5]/100;
dx = .0001;
N = ell/dx;
xgrid = [(1:N(1))'; (1:N(2))'; (1:N(3))']*dx;
NN = sum(N);
g_0 = 1/15;
g_L(1:N(1)) = linspace(0,1,N(1)).^2 + g_0;
g_L(N(1)+1:N(1)+N(2)) =  linspace(0,1,N(2)).^2 + g_0;
g_L(N(1)+N(2)+1:NN) =  (1-linspace(0,1,N(3))).^2 + g_0;

%form second difference matrix
S = diag(ones(NN-1,1),-1)+diag(ones(NN-1,1),1)+diag(-2*ones(NN,1));
S(1,1) = -1;
S(end) = -1;
S(N(1)+N(2)+1, N(1))=1;
S(N(1), N(1)+N(2)+1)=1;
S(N(1),N(1)+1) = 0;
S(N(1)+1,N(1)) = 0;
S(N(1)+1,N(1)+1)=-1;
S(N(1)+N(2)+1,N(1)+N(2)+1)=-3;
S = S*a(1)/dx/dx/2/R_a;

mat = S + diag(g_L);

[mvec mval] = eig(mat);

%store local q values
tq(1:numreq,1) = mvec(1,end-numreq+1:end)';
tq(1:numreq,2) = mvec(N(1)+1,end-numreq+1:end)';
tq(1:numreq,3) = mvec(end,end-numreq+1:end)';

%store eigenvalues
teigs = diag(mval(end-numreq+1:end,end-numreq+1:end));

subplot 131
plot(xgrid, g_L);
axis([0 max(xgrid) 0 2*max(g_L)])
title('g_L', 'fontsize', 15)

subplot 132
plot(mvec(:,end-3:end))
legend('4', '3', '2', '1')
evals = diag(mval(end-3:end,end-3:end));
title(num2str(evals))

subplot 133
[svec sval] = eig(S);
scaled = mvec(:,end-3:end)./svec(:,end-3:end);
plot(scaled)
legend('4', '3', '2', '1')
axis([0 NN -2 2])
scaled_eigval = evals ./ [1; (1:3)'.^2*pi^2/sum(ell)/sum(ell)];
title(num2str(scaled_eigval))

