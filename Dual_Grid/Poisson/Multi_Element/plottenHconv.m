clear all
close all
clc


i=0;
for H=[1 3 5 7]
i=i+1;
marker = 'o^svd+*<>';
c=0.3;

filename = ['Pconv_H' num2str(H) '_c' num2str(c) '.mat'];
load(filename);

semilogy(NrCellRange,errorL2,['-' marker(i) 'k'],'markerface','b')
hold on
% semilogy(NrCellRange,errorL2_interp,'--r','markerface','r')
% hold on
% semilogy(NrCellRange,abs(errorL2-errorL2_interp),'g')
xlim([NrCellRange(1) 20])
xlabel('N')
ylabel('L^2-error')
% title('potential \phi')
legend('interpolation','numerical',0)

clearallbut i

end

ylim([1e-10 10])
set(gca,'ytick',10.^(-12:2:2))






% i=0;
% for N=[1 3 5 7]
% i=i+1;
% marker = 'o^sv';
% c=0.3;
% 
% filename = ['Hconv_N' num2str(N) '_c' num2str(c) '.mat'];
% load(filename);
% 
% 
% HconvRange = 2.^(1:log2(Hconv));
% loglog(2./(HconvRange(2:end)),errorL2(2:end),['-' marker(i) 'b'],'markerface','b')
% hold on
% loglog(2./(HconvRange(2:end)),errorL2_interp(2:end),'--r')%,'markerface','r')
% hold on
% xlabel('h')
% ylabel('L^2-error \phi')
% % title('potential ')
% legend('interpolation','numerical',0)
% 
% ylim([1e-10 10])
% set(gca,'ytick',10.^(-12:2:2))
% 
% clearallbut i
% 
% end





