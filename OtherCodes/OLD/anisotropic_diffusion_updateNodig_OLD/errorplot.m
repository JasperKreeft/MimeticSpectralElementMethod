% figure
% semilogy(Z,errorL1,color(m))
% hold on
% xlim([Z(1) z])
% xlabel('N')
% ylabel('L^1-error')

figure
semilogy(Z,errorL2(Z))%,color(m))
hold on
semilogy(Z,errorL2_interp(Z),'--r')%['--' color(m)])
semilogy(Z,errorL2_uv(Z),'g')
xlim([Z(1) N])
xlabel('N')
ylabel('L^2-error')

legend('numerical','interpolation','velocity')