function [ts_g_L  G] = apply_q10(g_L, Temp, x_grid)
%function [ts_g_L  G] = apply_q10(g_L, Temp, x_grid)
%apply_q10 effectualizes temperature dependence of axial and
%membrane conductance.  
%
%input -
%  g_L : membrane leakage conductance at 35C
%  Temp : temperatures to find for
%
%output - 
%  ts_g_L : matrix of leakage conductances at temperatures 
%  G : matrix of axial conductance at temperatures (see Fohlmeister
%      et. al)


%find temperature-specific G_L from baseline g_L
q10_R_m = 1.85;
q10_g_L = 1.85;%  1/q10_R_m;
ref_temp = 35; %C

%use inverse q10 to get conductance
Qfacs_g_L = 1./(q10_g_L.^((ref_temp-Temp)/10));
ts_g_L = g_L(:)*(Qfacs_g_L(:)');


%find G, cytoplasmic conductivity
q10_R_i = .8;
ref_R_i = 140e-3; %Kohms*cm
ref_G = 1/ref_R_i;
ref_temp_i = 36; %C

q10_G = 1/q10_R_i;
Qfacs_G = 1./(q10_G.^((ref_temp_i-Temp)/10));
G = ref_G*Qfacs_G(:); %mS*cm

if nargin > 2
    if numel(x_grid)<2
        x_grid = 1:numel(g_L);
    end
    subplot 121
    [g1 g2] = meshgrid(Temp, x_grid);
    plot3(g1, g2, ts_g_L)
    xlabel('Temperature (C)', 'FontSize', 14)
    ylabel('Space (cm)', 'FontSize', 14)
    zlabel('mS/cm^2')
    title('Leakage Conductance', 'FontSize', 17)    

    subplot 122
    plot(Temp, G)
    xlabel('Temperature', 'FontSize', 14)
    ylabel('mS*cm', 'FontSize', 14)
    title('Axial Conductance', 'FontSize', 17)
end