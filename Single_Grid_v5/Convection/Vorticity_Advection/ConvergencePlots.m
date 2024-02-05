% close all
clc

format long e


kleur = 'brgc';
symbol = 'o^sv'
for N=1:4

    clearallbut kleur symbol N
    
    load(['LDC_Hconv_N' num2str(N) '_c2.mat'])
    
    ax = 1./(HconvRange);
    axisXY1 = [ 0.01 1 1e-10 1];
    axisXY2 = [ 0.01 1 1e-12 1];
%     str2 = kleur(N);
%     str3 = ['-' symbol(N) str2];
    str2 = 'w';
    str3 = ['--' symbol(N) kleur(N)];
    figure(11)
    loglog(ax,L2_U,str3,'markerface',str2)
    hold on
    axis(axisXY1)
    ylabel('||u-u_e_x||_{L^2\Lambda^1}')
    xlabel('h')
    set(gca,'ytick',10.^(-10:2:0))
    figure(12)
    loglog(ax,L2_w,str3,'markerface',str2)
    hold on
    axis(axisXY1)
    ylabel('||\omega-\omega_e_x||_{L^2\Lambda^0}')
    xlabel('h')
    set(gca,'ytick',10.^(-10:2:0))
    figure(13)
    loglog(ax,abs(Energy-Energy_exact),str3,'markerface',str2)
    hold on
    ylabel('\int_\Omega 1/2|u^2-u^2_{ex}|d\Omega')
    xlabel('h')
    axis(axisXY2)
    set(gca,'ytick',10.^(-12:2:0))
    figure(14)
    loglog(ax,abs(Enstrophy-Enstrophy_exact),str3,'markerface',str2)
    hold on
    ylabel('\int_\Omega 1/2|\omega^2-\omega^2_{ex}|d\Omega')
    xlabel('h')
    axis(axisXY2)
    set(gca,'ytick',10.^(-12:2:0))

end

format short
