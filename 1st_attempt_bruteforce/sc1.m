function [outG err] =  sc1 
close all 
cab = struct('Ri',0.06,'C',1,'a',1e-3,'L',0.1,'nocg',16,... 
             'nomult',4,'N',5000,'dt',2e-3); 
noc = cab.nocg*cab.nomult; 
dx = cab.L/(noc-1); 
G = (1:cab.nocg)'/cab.nocg; 
G = exp(-5*(G-.5).^2); 
%G = exp(G); 
%G = sqrt(G); 
%G = 0.5*ones(size(G)); 
for j=1:cab.nocg, 
    Ge(1+(j-1)*cab.nomult:j*cab.nomult) = G(j); 
end 
dB = -2*ones(noc,1); dB(1) = -1; dB(noc) = -1; 
oB = ones(noc,1); 
S = spdiags([oB dB oB],-1:1,noc,noc)*cab.a/(2*dx^2); 
t = 0:cab.dt:(cab.N-1)*cab.dt; 
i_0 = (sign(t-1) - sign(t-5))/800; 
i_0 = 5*t.*exp(-t)/800; 
V = zeros(noc,cab.N); 
VMEAS = zeros(5,cab.N); 
for th=1:5, 
    Rith = cab.Ri*th; 
    B = cab.dt*( S/Rith - spdiags(Ge',0,noc,noc)*th/2 )/cab.C; 
    SS = speye(noc) - B; 
    [LB,UB] = lu(SS); 
    for j=2:cab.N, 
        b = V(:,j-1); 
        b(1) = b(1) + cab.dt*i_0(j)/cab.C/(2*pi*cab.a*dx); 
        V(:,j) = UB \ ( LB \ b); 
    end 
    VMEAS(th,:) = V(1,:).*(1+randn(1,cab.N)/4000); 
end 
plot(t,VMEAS) 
drawnow 
legend('1','2','3','4','5') 
xlabel('t  (ms)','fontsize',14) 
ylabel('Soma potential (mV)','fontsize',14) 
options = optimset('display','iter',... 
                   'maxiter',100,... 
                   'gradobj','on',... 
                   'derivativecheck','off',... 
                   'largescale','off'); 
G_0 = mean(G)*ones(cab.nocg,1)+randn(cab.nocg,1)/50;

tic 
rG = fminunc(@(G) funth(G,cab,VMEAS),G_0,options); 
toc 
figure(2) 
x = linspace(0,cab.L,cab.nocg); 
plot(x,G); hold on; plot(x,rG,'r--'); plot(x,G_0,'-.o'); 
legend('true','recovered','initial') 
xlabel('x  (cm)','fontsize',14) 
ylabel('G_L (mS/cm^2)','fontsize',14) 
% 
% val and grad of cost function 
% 

for n = 1:cab.nocg
    outG((n-1)*cab.nomult+1::n*cab.nomult) = rG(n);
end
err = sum((outG(:) - Ge(:)).^2);


function [val,gval] = funth(G,cab,VMEAS) 
noc = cab.nocg*cab.nomult; 
dx = cab.L/(noc-1); 
gval = zeros(1,cab.nocg); 
for j=1:cab.nocg, 
    Ge(1+(j-1)*cab.nomult:j*cab.nomult) = G(j); 
end 
% figure(2) 
% plot(G,'rx--') 
% drawnow 
dB = -2*ones(noc,1); dB(1) = -1; dB(noc) = -1; 
oB = ones(noc,1); 
S = spdiags([oB dB oB],-1:1,noc,noc); 
t = 0:cab.dt:cab.N*cab.dt; 
%i_0 = (sign(t-1) - sign(t-5))/800/C/(2*pi*a*dx); 
i_0 = 5*t.*exp(-t)/800/cab.C/(2*pi*cab.a*dx); 
val = 0; 
for th=1:5, 
    v = zeros(noc,cab.N); 
    B = cab.dt*( (cab.a/(2*cab.Ri*th*dx^2))*S - ... 
              spdiags(Ge',0,noc,noc)*th/2 )/cab.C; 
    SS = speye(noc) - B; 
    [LB,UB] = lu(SS); 
    for j=2:cab.N, 
        b = v(:,j-1); 
        b(1) = b(1) + cab.dt*i_0(j); %/C/(2*pi*a*dx); 
        v(:,j) = UB \ ( LB \ b); 
    end 
    val = val + sum((v(1,:)-VMEAS(th,:)).^2)*cab.dt/2; 
    % solve the adjoint problem 
    V = zeros(noc,cab.N); 
    for j=cab.N-1:-1:1 
        b = V(:,j+1); 
        b(1) = b(1) + cab.dt*(VMEAS(th,j)-v(1,j))/dx; 
        V(:,j) = UB \ ( LB \ b);

    end 
    % evaluate the gradient 
    for i=1:cab.nocg, 
        xind = 1+(i-1)*cab.nomult:i*cab.nomult; 
        gval(i) = gval(i) + ... 
                  sum(sum(V(xind,:).*v(xind,:)))*cab.dt*dx*(th/2); 
    end 
end 
