clear all
close all
clc

errorplot = false;

if errorplot == true
    error_L2 = []; error_L2_dudx = []; hh = [];
    Ne = [4 6 8 10 15 20 25];         % Nr of elements
else
    Ne = 8;         % Nr of elements
end

for Ne = Ne

    Nn = Ne+1;      % Nr of nodes
    L  = 1;         % length domain

    % equidistant space
    % x = linspace(0,L,Ne+1);
    % one-sided refinement
    ratio = 1/1;
    x = elm_line1 (0,L,Ne,ratio);
    h = diff(x);

    bc_R = 0;
    bc_L = 0;

    % % Exact solution
    if errorplot == false
        figure(1)
        xx = linspace(0,1,100);
        yy = (bc_L*exp(1)-bc_R+1)/(exp(1)-exp(-1))*exp(-xx)+(bc_R-bc_L*exp(-1)-1)/(exp(1)-exp(-1))*exp(xx)+xx;
        plot(xx,yy,'g','linewidth',2)
        figure(2)
        yy = -(bc_L*exp(1)-bc_R+1)/(exp(1)-exp(-1))*exp(-xx)+(bc_R-bc_L*exp(-1)-1)/(exp(1)-exp(-1))*exp(xx)+ones(1,length(xx));
        plot(xx,yy,'g','linewidth',2)
    end

    GM = zeros(Ne,2);
    GM(:,1) = (1:Ne)';
    GM(:,2) = (1:Ne)'+1;

    F = zeros(Nn,1);
    K = zeros(Nn);
    K_u = zeros(Nn-1,1);
    K_d = zeros(Nn,1);
    K_l = zeros(Nn-1,1);

    for k=1:Ne

        Kk11 = 1/h(k)+h(k)/3; Kk12 = -1/h(k)+h(k)/6;
        Kk21 = Kk12;          Kk22 = 1/h(k)+h(k)/3;

    %     K_u(k)   = Kk12;
    %     K_d(k)   = K_d(k)+Kk11;
    %     K_d(k+1) = K_d(k+1)+Kk22;
    %     K_l(k)   = Kk21;

        K(k:k+1,k:k+1) = K(k:k+1,k:k+1)+[Kk11 Kk12; Kk21 Kk22];

        Fk = h(k)/6*[2*x(k)+x(k+1); x(k)+2*x(k+1)];
        for i=1:2
            F(GM(k,i)) = F(GM(k,i))+Fk(i);
        end
    end

    F(2) = F(2)-bc_L*K(2,1);
    F(end-1) = F(end-1)-bc_R*K(end-1,end);

    A = K(2:end-1,2:end-1)\F(2:end-1);
    % A = thomas (Nn-2,K_d(2:end-1),K_u(2:end),K_l(1:end-1),F(2:end-1));
    A = [bc_L;A;bc_R];

    if errorplot == false
        figure(1)
        hold on
        plot(x,A,'-o','linewidth',2)
        grid on
        set(gca,'FontSize',14)
        xlabel('x','FontSize',16)
        ylabel('u','FontSize',16)
        legend('exact','FEM',2)
        xlim([0 1])


%         YI = interp1(x,A',linspace(x(1),x(end),100),'spline');
%         plot(linspace(x(1),x(end),100),YI,'r')

        figure(2)
        hold on
        dAdx = diff(A)'./h;
        stairs(x,[dAdx dAdx(end)],'linewidth',2);
        grid on
        set(gca,'FontSize',14)
        xlabel('x','FontSize',16)
        ylabel('du/dx','FontSize',16)
        legend('exact','FEM',3)
        axis([0 1 -0.4 0.2])

    else

        if roundn(h-h(1),-14)==0
            xxx = linspace(0,1,1000*Ne+1);
            for k=1:Ne
                i=k+1;
                for j=(1:1001)+(k-1)*1000
                    u_app(j)    = A(i-1)+(xxx(j)-x(i-1))/h(k)*(A(i)-A(i-1));
                    dudx_app(j) = (A(i)-A(i-1))/h(k);
                    u_ex(j)     = (bc_L*exp(1)-bc_R+1)/(exp(1)-exp(-1))*exp(-xxx(j))+(bc_R-bc_L*exp(-1)-1)/(exp(1)-exp(-1))*exp(xxx(j))+xxx(j);
                    dudx_ex(j)  = -(bc_L*exp(1)-bc_R+1)/(exp(1)-exp(-1))*exp(-xxx(j))+(bc_R-bc_L*exp(-1)-1)/(exp(1)-exp(-1))*exp(xxx(j))+1;
                end
            end
            error_L2 = [error_L2 ; sqrt(sum((u_ex-u_app).^2*(xxx(2)-xxx(1))))];
            error_L2_dudx = [error_L2_dudx ; sqrt(sum((dudx_ex-dudx_app).^2*(xxx(2)-xxx(1))))];
            hh = [hh ; 1/Ne];
        end

    end

end

if errorplot == true

    figure(3)
    loglog(hh,error_L2,'-o','linewidth',2)
    line([7/9*1e-1 1e-1 1e-1 7/9*1e-1],[exp(-8.5172) exp(-7.9852) exp(-8.5172) exp(-8.5172)],'color','k')
    axis([3e-2 3e-1 6e-5 4e-3])
    grid on
    set(gca,'FontSize',14)
    xlabel('h','FontSize',16)
    ylabel('e','FontSize',16)
    text(9e-2,exp(-8.4),'2','FontSize',14)

    figure(4)
    loglog(hh,error_L2_dudx,'-o','linewidth',2)
    line([7/9*1e-1 1e-1 1e-1 7/9*1e-1],[exp(-4.6052) exp(-4.3275) exp(-4.6052) exp(-4.6052)],'color','k')
    axis([3e-2 3e-1 5e-3 5e-2])
    grid on
    set(gca,'FontSize',14)
    xlabel('h','FontSize',16)
    ylabel('de/dx','FontSize',16)
    text(9e-2,exp(-4.55),'1','FontSize',14)

end