%    stim_int(:,time<=t1) = 0;
    stim_on = (time>=t1).*(time<t2);
%    ib_time = time(stim_on);
%    repzeta = repmat(sc_zeta,Nx,numel(ib_time));
    
%    stim_int(:,stim_on) = (exp(repzeta.*(ib_time-t1))-1)./repzeta;


%first rearrangement
    %third check: rearrangement
    top_l = V(n,1).*(Gprime(Tgrid(1))/dx/dx.*(V(n+1,1))...
        + G(Tgrid(1))/dx/dx/dtemp.*(V(n+1,2)-V(n+1,1)))...
        - V(n+1,1)*G(Tgrid(1))/dx/dx/dtemp.*(V(n,2)-V(n,1));
    
    mid_l = V(n,2:end-1).*(Gprime(itv)/dx/dx.*(V(n+1,2:end-1))...
        + G(itv)/dx/dx/2/dtemp.*(V(n+1,3:end)-V(n+1,1:end-2))) ...
        - V(n+1,2:end-1).*G(itv)/dx/dx/2/dtemp.*(V(n,3:end)-V(n,1:end-2)) ;
    
    bot_l = V(n,end).*(Gprime(Tgrid(end))/dx/dx.*(V(n+1,end))...
        + G(Tgrid(end))/dx/dx/dtemp.*(V(n+1,end)-V(n+1,end-1)))...
        - V(n+1,end).*G(Tgrid(end))/dx/dx/dtemp.*(V(n,end)-V(n,end-1));
    
    check_lhs3 = [top_l mid_l bot_l];
    
    
    top_r = G(Tgrid(1))/dx/dx/dtemp.*(-2*V(n,1)+V(n-1,1)).*(V(n,2)-V(n,1))...
        -V(n,1).*(G(Tgrid(1))/dx/dx/dtemp.*(-2*V(n,2)+V(n-1,2)+2*V(n,1)-V(n-1,1))...
                + Gprime(Tgrid(1))/dx/dx.*(-2*V(n,1)+V(n-1,1)));
    
    mid_r = G(itv)/dx/dx/2/dtemp.*(-2*V(n,2:end-1)+V(n-1,2:end-1)).*(V(n,3:end)-V(n,1:end-2))...
        -V(n,2:end-1).*(G(itv)/dx/dx/2/dtemp.*(-2*V(n,3:end)+V(n-1,3:end)+2*V(n,1:end-2)-V(n-1,1:end-2))...
                      + Gprime(itv)/dx/dx.*(-2*V(n,2:end-1)+V(n-1,2:end-1)));
    
    bot_r = G(Tgrid(end))/dx/dx/dtemp.*(-2*V(n,end)+V(n-1,end)).*(V(n,end)-V(n,end-1))...
        -V(n,end).*(G(Tgrid(end))/dx/dx/dtemp.*(-2*V(n,end)+V(n-1,end)+2*V(n,end-1)-V(n-1,end-1))...
                  + Gprime(Tgrid(end))/dx/dx.*(-2*V(n,end)+V(n-1,end)));

    check_rhs3 = [top_r mid_r bot_r];
    
    
    
%    more checks:


    %first, check our already known solution
    vxx = (V(n+1,:)-2*V(n,:)+V(n-1,:))/dx/dx;
    vxxt(2:Ntemp-1) = (vxx(3:end)-vxx(1:end-2))/2/dtemp;
    vxxt(1) = (vxx(2)-vxx(1))/dtemp;
    vxxt(Ntemp) = (vxx(Ntemp)-vxx(Ntemp-1))/dtemp;
    vt(2:Ntemp-1) = (V(n,3:Ntemp)-V(n,1:Ntemp-2))/2/dtemp;
    vt(1) = (V(n,2)-V(n,1))/dtemp;
    vt(Ntemp) = (V(n,Ntemp)-V(n,Ntemp-1))/dtemp;

    check_lhs = V(n,:).*(Gprime(Tgrid).*vxx + G(Tgrid).*vxxt);
    check_rhs = G(Tgrid).*vxx.*vt;

    
    %second check
    top_l = V(n,1).*(Gprime(Tgrid(1))/dx/dx.*(V(n+1,1)-2*V(n,1)+V(n-1,1))...
        + G(Tgrid(1))/dx/dx/dtemp.*(V(n+1,2)-2*V(n,2)+V(n-1,2)-V(n+1,1)+2*V(n,1)-V(n-1,1)));
    mid_l = V(n,2:end-1).*(Gprime(itv)/dx/dx.*(V(n+1,2:end-1)-2*V(n,2:end-1)+V(n-1,2:end-1))...
        + G(itv)/dx/dx/2/dtemp.*(V(n+1,3:end)-2*V(n,3:end)+V(n-1,3:end)-V(n+1,1:end-2)+2*V(n,1:end-2)-V(n-1,1:end-2)));
    bot_l = V(n,end).*(Gprime(Tgrid(end))/dx/dx.*(V(n+1,end)-2*V(n,end)+V(n-1,end))...
        + G(Tgrid(end))/dx/dx/dtemp.*(V(n+1,end)-2*V(n,end)+V(n-1,end)-V(n+1,end-1)+2*V(n,end-1)-V(n-1,end-1)));
    
    check_lhs2 = [top_l mid_l bot_l];
    
    top_r = G(Tgrid(1))/dx/dx/dtemp.*(V(n+1,1)-2*V(n,1)+V(n-1,1)).*(V(n,2)-V(n,1));
    mid_r = G(itv)/dx/dx/2/dtemp.*(V(n+1,2:end-1)-2*V(n,2:end-1)+V(n-1,2:end-1)).*(V(n,3:end)-V(n,1:end-2));
    bot_r = G(Tgrid(end))/dx/dx/dtemp.*(V(n+1,end)-2*V(n,end)+V(n-1,end)).*(V(n,end)-V(n,end-1));

    check_rhs2 = [top_r mid_r bot_r];
    
    
    
    
    %third check: rearrangement
    top_l = V(n+1,1).*(V(n,1).*(Gprime(Tgrid(1))/dx/dx - G(Tgrid(1))/dx/dx/dtemp) - G(Tgrid(1))/dx/dx/dtemp.*(V(n,2)-V(n,1)))...
        + V(n+1,2).*(V(n,1).*(G(Tgrid(1))/dx/dx/dtemp));
    
    mid_l = V(n+1,2:end-1).*(V(n,2:end-1).*(Gprime(itv)/dx/dx) - G(itv)/dx/dx/2/dtemp.*(V(n,3:end)-V(n,1:end-2)))...
        + V(n+1,3:end).*(V(n,2:end-1).*(G(itv)/dx/dx/2/dtemp))...
        +V(n+1,1:end-2).*(-V(n,2:end-1).*(G(itv)/dx/dx/2/dtemp));
    
    bot_l = V(n+1,end).*(V(n,end).*(Gprime(Tgrid(end))/dx/dx + G(Tgrid(end))/dx/dx/dtemp )-G(Tgrid(end))/dx/dx/dtemp.*(V(n,end)-V(n,end-1)))...
        + V(n+1,end-1)*(-V(n,end).*(G(Tgrid(end))/dx/dx/dtemp));
    
    
    check_lhs3 = [top_l mid_l bot_l];
    
    
    top_r = G(Tgrid(1))/dx/dx/dtemp.*(-2*V(n,1)+V(n-1,1)).*(V(n,2)-V(n,1))...
        -V(n,1).*(G(Tgrid(1))/dx/dx/dtemp.*(-2*V(n,2)+V(n-1,2)+2*V(n,1)-V(n-1,1))...
                + Gprime(Tgrid(1))/dx/dx.*(-2*V(n,1)+V(n-1,1)));
    
    mid_r = G(itv)/dx/dx/2/dtemp.*(-2*V(n,2:end-1)+V(n-1,2:end-1)).*(V(n,3:end)-V(n,1:end-2))...
        -V(n,2:end-1).*(G(itv)/dx/dx/2/dtemp.*(-2*V(n,3:end)+V(n-1,3:end)+2*V(n,1:end-2)-V(n-1,1:end-2))...
                      + Gprime(itv)/dx/dx.*(-2*V(n,2:end-1)+V(n-1,2:end-1)));
    
    bot_r = G(Tgrid(end))/dx/dx/dtemp.*(-2*V(n,end)+V(n-1,end)).*(V(n,end)-V(n,end-1))...
        -V(n,end).*(G(Tgrid(end))/dx/dx/dtemp.*(-2*V(n,end)+V(n-1,end)+2*V(n,end-1)-V(n-1,end-1))...
                  + Gprime(Tgrid(end))/dx/dx.*(-2*V(n,end)+V(n-1,end)));

    check_rhs3 = [top_r mid_r bot_r]
    
    
    
    
%    the buggy right-hand side:

    rhs(1) = sV(n,1)*(Gprime(mintemp)/dx/dx*(2*sV(n,1) + sV(n-1,1))...
        + G(mintemp)/dtemp/dx/dx*(2*sV(n,2)-sV(n-1,2)-2*sV(n,1)+sV(n-1,1) ) )...
        + G(mintemp)/dtemp/dx/dx*(-2*sV(n,1)+sV(n-1,1))*(sV(n,2)-sV(n,1));

    rhs(2:Ntemp-1) = (sV(n,2:end-1)/dx/dx.*Gprime(itv).*(2*sV(n,2:end-1)-sV(n-1,2:end-1))...
        + G(itv)/2/dx/dx/dtemp.*(2*sV(n,3:end)-sV(n-1,3:end)-2*sV(n,1:end-2)+sV(n-1,1:end-2)) ...
        +G(itv)/2/dx/dx/dtemp.*(-2*sV(n,2:end-1)+sV(n-1,2:end-1)).*(sV(n,3:end)-sV(n,1:end-2)))';

    rhs(Ntemp) = sV(n,end)*(Gprime(Tgrid(Ntemp))/dx/dx*(2*sV(n,end) + sV(n-1,end))...
        + G(mintemp)/dtemp/dx/dx*(2*sV(n,end)-sV(n-1,end)-2*sV(n,end-1)+sV(n-1,end-1) ) )...
        + G(mintemp)/dtemp/dx/dx*(-2*sV(n,end)+sV(n-1,end))*(sV(n,end)-sV(n,end-1));