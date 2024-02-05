% Curves problem
clear all
close all
hold on
clc

% k1 = clock;

T  = 2*pi;
dt0 = pi;


disp('Duffing equation problem')
disp('  ')
disp(' 1. Backward Euler')
disp(' 2. BDF 2nd order')
disp(' 3. Trapeziodal method')
disp(' 4. ESDIRK 3rd order')
disp(' 5. ESDIRK 4th order')
disp(' 6. ESDIRK 5th order')
disp('  ')
methods = sort(str2num(input('Enter schemes to be displayed (separated by spaces): ','s')));

legend_strings = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Backward Euler
if methods(1) == 1

    methods(1) = [];
    if isempty(methods)
        methods = 10;
    end

    dt = dt0/7;
    dtau = dt;
    N    = floor(T/dt)+1;

    t    = zeros(1,N);
    X_BE = zeros(2,N);

    X_BE(:,1) = [-0.2;0.2];

    for n=1:N-1
        t(n+1) = t(n)+dt;
        Xn_new = X_BE(:,n);
        Xn_old = [10;10];
        while abs(max(Xn_new-Xn_old)./Xn_old)>1e-12
            Xn_old = Xn_new;
            [A] = vector_elements(Xn_old);
            Xn_new = Xn_old-dtau*((Xn_old-X_BE(:,n))/dt-A*Xn_old);
        end
        X_BE(:,n+1) = Xn_new;
    end
    plot(X_BE(1,:),X_BE(2,:),'.-k','linewidth',2)
    legend_strings = [legend_strings;'BE     '];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BDF2
if methods(1) == 2

    methods(1) = [];
    if isempty(methods)
        methods = 10;
    end

    dt = dt0/7;
    dtau = dt;
    N    = floor(T/dt)+1;

    t      = zeros(1,N);
    X_BDF2 = zeros(2,N);

    X_BDF2(:,1) = [-0.2;0.2];

    for n=1:N-1
        t(n+1) = t(n)+dt;

        i=0;
        Xn_new = X_BDF2(:,n);
        Xn_old = [10;10];
        while abs(max(Xn_new-Xn_old))>1e-12
            Xn_old = Xn_new;
            [A] = vector_elements(Xn_old);
            if n==1
                Xn_new = Xn_old-dtau*((Xn_old-X_BDF2(:,n))/dt-A*Xn_old);
            else
                Xn_new = Xn_old-dtau*((3/2*Xn_old-2*X_BDF2(:,n)+1/2*X_BDF2(:,n-1))/dt-A*Xn_old);
            end
        end
        X_BDF2(:,n+1) = Xn_new;
    end
    plot(X_BDF2(1,:),X_BDF2(2,:),'.-m','linewidth',2)
    legend_strings = [legend_strings;'BDF2   '];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if methods(1) == 3
    %%% CN / trapeziodal with pseudo time stepping

    methods(1) = [];
    if isempty(methods)
        methods = 10;
    end

    dt = dt0/7;
    dtau = dt;
    N    = floor(T/dt)+1;

    t    = zeros(1,N);
    X_CN = zeros(2,N);

    X_CN(:,1) = [-0.2 ; 0.2];

    for n=1:N-1
        t(n+1) = t(n)+dt;
        [A_n] = vector_elements(X_CN(:,n));

        Xn_new = X_CN(:,n);
        Xn_old = [10;10];
        while abs(max(Xn_new-Xn_old)./Xn_old)>1e-12
            Xn_old = Xn_new;
            [A_np1] = vector_elements(Xn_old);
            Xn_new = Xn_old-dtau*((Xn_old-X_CN(:,n))/dt-1/2*(A_np1*Xn_old)-1/2*(A_n*X_CN(:,n)));
        end
        X_CN(:,n+1) = Xn_new;
    end
    legend_strings = [legend_strings;'Trap/CN'];
    plot(X_CN(1,:),X_CN(2,:),'.-c','linewidth',2)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while methods(1)<=6

if methods(1) == 4 || methods(1) == 5 || methods(1) == 6

    if methods(1) == 4
        dt = dt0*3/7;
        esdirk3;
        kleur = 'g';
        legend_strings = [legend_strings;'ESDIRK3'];
    elseif methods(1) == 5
        dt = dt0*5/7;
        esdirk4;
        kleur = 'b';
        legend_strings = [legend_strings;'ESDIRK4'];
    elseif methods(1) == 6
        dt = dt0;
        esdirk5;
        kleur = 'r';
        legend_strings = [legend_strings;'ESDIRK5'];
    end

    methods(1) = [];
    if isempty(methods)
        methods = 10;
    end

    dtau = dt;
    N        = floor(T/dt)+1;
    t        = zeros(1,N);
    X_ESDIRK = zeros(2,N);

    X_ESDIRK(:,1) = [-0.2;0.2];

    s = length(rkaI);

    for n=1:N-1

        Xe = zeros(2,s);
        Xn_new = zeros(2,s);
        Xn_old = zeros(2,s);


        t(n+1) = t(n)+dt;

        k=1;
        Xe(:,1) = X_ESDIRK(:,n);

        for k=2:s

            S = [0;0];
            for j=1:k-1
                [A_j] = vector_elements(Xe(:,j));
                S = S+rkaI(k,j)*(A_j*Xe(:,j));
            end

            Xn_new = X_ESDIRK(:,n);
            Xn_old = [10;10];
            err = 1;
            while abs(max((Xn_new-Xn_old)./Xn_old))>1e-12
                Xn_old = Xn_new;
                [A_k] = vector_elements(Xn_old);
                Xn_new = Xn_old-dtau*((Xn_old-X_ESDIRK(:,n))/dt-rkaI(k,k)*(A_k*Xn_old)-S);
    % SvZ       Xn_new = (eye(2)-dt*rka(k,k)*A_k)\(X_ESDIRK(:,n)+dt*S);
            end
            Xe(:,k) = Xn_new;
        end
        X_ESDIRK(:,n+1) = Xe(:,s);
    end
    plot(X_ESDIRK(1,:),X_ESDIRK(2,:),['.-' kleur],'linewidth',2)
end

end %while
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

legend(legend_strings,0)
grid on
title('Curves')
ylabel('X_1')
xlabel('time t')