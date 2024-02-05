load divergence_free.mat

% N=2, Cartesian mesh.

figure
loglog(ax(2:end),L2_du(2:end),'^k','markerface','k')
hold on
loglog(ax(2:end),L1_du(2:end),'vk','markerface','k')
loglog(ax(2:end),Linf_du(2:end),'ok','markerface','k')
axis([0.02 2 1e-16 1e-13])
title('pointwise divergence-free solution')
ylabel('L^1, L^2, L^\infty-error')
xlabel('h')
legend('L^2-error','L^1-error','L^\infty-error',1)