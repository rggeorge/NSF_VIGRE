if 1
    Ntemp = [1 1 1 1 1];

    clear rG err G

    [rG(:,1) err(1) G f_err(1)] = sc2(Ntemp(1));

    for n = 2:numel(Ntemp)
        [rG(:,n) err(n) G f_err(n)] = sc2(Ntemp(n),rG(:,n-1) );
    end
end

figure(5)
subplot 121
plot(Ntemp, err);

subplot 122
plot(rG)
hold on
plot(G, 'Linewidth', 2)
legend('1', '2', '3', '4', '5', '6', '7', 'True g_L')
