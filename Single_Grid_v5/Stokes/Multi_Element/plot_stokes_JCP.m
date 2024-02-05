close all

% kleur = 'brgc';
% mark = 'o^sv';
% axisXY = [0.009 1 1e-12 1];
% 
% i=0;
% for N=[2 4 6 8]
%     i=i+1;
%     
%     load(['LDC_Hconv_N' num2str(N) '_c2.mat'])
%     
%     ax = 1./HconvRange;
%     
%     str1 = ['-' mark(i) kleur(i)];
%     str2 = kleur(i);
%     str3 = ['--' kleur(i)];
%     
%     figure(1)
%     loglog(ax,Hd_w,str1,'markerface',str2)
%     hold on
%     loglog(ax,L2_w,str3,'markerface',str2)
%     axis(axisXY);
%     set(gca,'ytick',10.^(-12:2:0))
%     
%     figure(2)
%     loglog(ax,Hd_u,str1,'markerface',str2)
%     hold on
%     loglog(ax,L2_U,str3,'markerface',str2)
%     axis(axisXY);
%     set(gca,'ytick',10.^(-12:2:0))
%     
%     figure(3)
%     loglog(ax,L2_p,str1,'markerface',str2)
%     hold on
%     axis(axisXY);
%     set(gca,'ytick',10.^(-12:2:0))
%     
%     
% end
% 
% figure(4)
% for r=1:9
% Y = exp(r*log([1/2 1/64]));
% loglog([1/2 1/64],Y,'k')
% hold on
% end
% axis(axisXY)
% set(gca,'ytick',10.^(-12:2:0))


kleur = 'brgc';
mark = 'o^sv';
axisXY = [ 1 17 1e-12 1];

i=0;
for Hconv=[1 2 4 8]
    i=i+1;q

    load(['LDC_Pconv_H' num2str(Hconv) '_c0.mat'])

    ax = NrCellRange;

    str1 = ['-' mark(i) kleur(i)];
    str2 = kleur(i);
    str3 = ['--' kleur(i)];

    figure(1)
    semilogy(ax,Hd_w,str1,'markerface',str2)
    hold on
    semilogy(ax,L2_w,str3,'markerface',str2)
    axis(axisXY)
    set(gca,'ytick',10.^(-12:2:0))
    set(gca,'xtick',2:2:16)
    figure(2)
    semilogy(ax,L2_u,str1,'markerface',str2)
    hold on
    axis(axisXY)
    set(gca,'ytick',10.^(-12:2:0))
    set(gca,'xtick',2:2:16)
    figure(3)
    semilogy(ax,L2_p,str1,'markerface',str2)
    hold on
    axis(axisXY)
    set(gca,'ytick',10.^(-12:2:0))
    set(gca,'xtick',2:2:16)

end
