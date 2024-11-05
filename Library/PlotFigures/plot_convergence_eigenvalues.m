function plot_convergence_eigenvalues(NH,range,error,exact)

global N

switch NH
    case "N"

if length(range)>=2
    figure
    handle(9) = Marker(semilogy(range,error(9,:),'-o'));%,':dk','markerface','k');
    hold on
    handle(8) = Marker(semilogy(range,error(8,:),'-o'));%,'--sk','markerface','k');
    handle(7) = Marker(semilogy(range,error(7,:),'-o'));%,'-ok','markerface','k');
    handle(6) = Marker(semilogy(range,error(6,:),'-o'));%,'-oy','markerface','y');
    handle(5) = Marker(semilogy(range,error(5,:),'-o'));%,'-oc','markerface','c');
    handle(4) = Marker(semilogy(range,error(4,:),'-o'));%,'-om','markerface','m');
    handle(3) = Marker(semilogy(range,error(3,:),'-o'));%,'-or','markerface','r');
    handle(2) = Marker(semilogy(range,error(2,:),'-o'));%,'-og','markerface','g');
    handle(1) = Marker(semilogy(range,error(1,:),'-o'));%,'-ob','markerface','b');
    grid on
    exact = string(exact(1:9));
    legend(handle,exact,'location','southwest')%,'orientation','horizontal')
    axis([0 N 1e-10 1e2])
    xlabel('N')
    ylabel('error eigenvalues')
    title('Convergence of first nine non-zero eigenvalues')
end


% if length(range)>=9
% figure(2)
% handle(9) = semilogy(range,error(9,:)',':dk','markerface','k');
% hold on
% handle(8) = semilogy(range,error(8,:)','--sk','markerface','k');
% handle(7) = semilogy(range,error(7,:)','-ok','markerface','k');
% handle(6) = semilogy(range,error(6,:)','-oy','markerface','y');
% handle(5) = semilogy(range,error(5,:)','-oc','markerface','c');
% handle(4) = semilogy(range,error(4,:)','-om','markerface','m');
% handle(3) = semilogy(range,error(3,:)','-or','markerface','r');
% handle(2) = semilogy(range,error(2,:)','-og','markerface','g');
% handle(1) = semilogy(range,error(1,:)','-ob','markerface','b');
% grid on
% legend(handle,'1','4','9','16','25','36','49','64','81','location','southwest')%,'orientation','horizontal')
% axis([0 N 1e-10 1e2])
% xlabel('N')
% ylabel('error eigenvalues')
% title('Convergence of first nine non-zero eigenvalues')
% end

    case "H"
   
        nrH = length(range);
        handle(9) = Marker(loglog(1./range,error(9,1:nrH),'-o'));
        hold on
        handle(8) = Marker(loglog(1./range,error(8,1:nrH),'-o'));
        handle(7) = Marker(loglog(1./range,error(7,1:nrH),'-o'));
        handle(6) = Marker(loglog(1./range,error(6,1:nrH),'-o'));
        handle(5) = Marker(loglog(1./range,error(5,1:nrH),'-o'));
        handle(4) = Marker(loglog(1./range,error(4,1:nrH),'-o'));
        handle(3) = Marker(loglog(1./range,error(3,1:nrH),'-o'));
        handle(2) = Marker(loglog(1./range,error(2,1:nrH),'-o'));
        handle(1) = Marker(loglog(1./range,error(1,1:nrH),'-o'));
        Conv = range.^(-2*N);
        C = 0.01;
        loglog(1./range(2:5),C*Conv(2:5),'k')
        loglog([1/range(2) 1/range(2)],[C*Conv(2) C*Conv(5)],'k')
        loglog([1/range(2) 1/range(5)],[C*Conv(5) C*Conv(5)],'k')
        Rate = (log(error(1,4))-log(error(1,2)))/(log(1/range(4))-log(1/range(2)));
        text(2/sum(range(3:4)),C*Conv(4),num2str(Rate,'% 2.1f'),'FontSize',18)
        % grid on
        legend(handle,'1','4','9','16','25','36','49','64','81','location','northwest')%,'orientation','horizontal')
        axis([1e-4 1 1e-10 1e2])
        xlabel('h')
        ylabel('error eigenvalues')
        title('h-Convergence of first nine non-zero eigenvalues')

end