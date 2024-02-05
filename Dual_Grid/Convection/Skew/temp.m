clear all
close all
clc

disp(' 1. Forward Euler')
disp(' 2. Backward Euler')
disp(' 3. Trapeziodal')
disp(' 4. BDF2')
disp(' 5. ERK3')
disp(' 6. ERK4')
disp(' 7. IRK2')
disp(' 8. ESDIRK3')
disp(' 9. ESDIRK4')
disp('10. ESDIRK5')

methods = sort(str2num(input('Enter schemes to be displayed (separated by spaces): ','s')));


T=60;
dt=2*pi/7;
N=round(T/dt);
t   = (0:N-1)*dt;

k  = 1;
m  = 1;
wn = sqrt(k/m);

R = [ 0 1 ; -wn^2 0 ];


% Exact solution

t_ex = 0:0.01:T;
q_ex = cos(wn*t_ex);
plot(t_ex,q_ex,'k')
hold on
axis([0 60 -2 2])
grid on
legend_strings = 'Exact  ';


if methods(1) == 1
    % Forward Euler 
    
    methods(1) = [];
    if isempty(methods)
        methods = 99;
    end
    
    Q   = zeros(2,N);
    
    Q(:,1) = [1;0];

    S = (eye(2)+dt*R);
    
    for n=2:N
        Q(:,n) = S*Q(:,n-1);
    end

    q_forw = Q(1,:);

    plot(t,q_forw,'.-g')

    legend_strings = [legend_strings;'FE     '];
end

if methods(1) == 2
    % Backward Euler

    methods(1) = [];
    if isempty(methods)
        methods = 99;
    end

    Q   = zeros(2,N);

    Q(:,1)   = [1;0];

    S = inv(eye(2)-dt*R);
    
    for n=2:N
        Q(:,n) = S*Q(:,n-1);
    end

    q_back = Q(1,:);

    plot(t,q_back,'.-r')

    legend_strings = [legend_strings;'BE     '];

end

if methods(1) == 3
    % Trapeziodal

    methods(1) = [];
    if isempty(methods)
        methods = 99;
    end

    Q   = zeros(2,N);

    Q(:,1)   = [1;0];
    
    theta = 1/2;

    S = inv(eye(2)-theta*dt*R)*(eye(2)+(1-theta)*dt*R);
    
    for n=2:N
        Q(:,n) = S*Q(:,n-1);
    end

    q_trap = Q(1,:);

    plot(t,q_trap,'.-')

    legend_strings = [legend_strings;'Trap/CN'];
end

if methods(1) == 4
    % BDF2

    methods(1) = [];
    if isempty(methods)
        methods = 99;
    end

    Q   = zeros(2,N);

    Q(:,1) = [1;0];

    Q(:,2)   = [cos(dt);-sin(dt)];

    S = inv(3/2*eye(2)-dt*R);
    
    for n=3:N
        Q(:,n) = S*(2*Q(:,n-1)-1/2*Q(:,n-2));
    end

    q_BDF2 = Q(1,:);

    plot(t,q_BDF2,'.-m')

    legend_strings = [legend_strings;'BDF2   '];
end

if methods(1) == 5
    % ERK3

    methods(1) = [];
    if isempty(methods)
        methods = 99;
    end
    
    Q   = zeros(2,N);

    Q(:,1) = [1;0];

    A = [  0   0   0 
           1   0   0 
          1/4 1/4  0];

    b = [1/6 1/6 2/3];

    n_a = length(A);
    
    AA = zeros(2*n_a);
    for i=1:2*n_a
        if mod(i,2)~=0
            AA(i,1:2:2*n_a) = A(ceil(i/2),:);
        else
            AA(i,2:2:2*n_a) = A(ceil(i/2),:);
        end
    end
    
    RR = zeros(2*n_a);
    for i=1:n_a
        RR(2*i-1:2*i,2*i-1:2*i) = R;
    end

    S = inv(inv(RR)-dt*AA);
    
    for n=2:N
        Qprime = S*[Q(:,n-1);Q(:,n-1);Q(:,n-1)];
        Q(1,n) = Q(1,n-1)+dt*b*Qprime(1:2:2*n_a);
        Q(2,n) = Q(2,n-1)+dt*b*Qprime(2:2:2*n_a);
    end

    q_ERK3 = Q(1,:);

    plot(t,q_ERK3,'.-c')

    legend_strings = [legend_strings;'ERK3   '];
end

if methods(1) == 6
    % ERK4

    methods(1) = [];
    if isempty(methods)
        methods = 99;
    end

    Q   = zeros(2,N);
    
    Q(:,1) = [1;0];
    
    A = [  0   0   0   0
          1/2  0   0   0
           0  1/2  0   0
           0   0   1  0];

    b = [1/6 1/3 1/3 1/6];

    n_a = length(A);
    
    AA = zeros(2*n_a);
    
    for i=1:2*n_a
        if mod(i,2)~=0
            AA(i,1:2:2*n_a) = A(ceil(i/2),:);
        else
            AA(i,2:2:2*n_a) = A(ceil(i/2),:);
        end
    end
    
    RR = zeros(2*n_a);
    for i=1:n_a
        RR(2*i-1:2*i,2*i-1:2*i) = R;
    end
    
    S = inv(inv(RR)-dt*AA);
    
    for n=2:N
        Qprime = S*[Q(:,n-1);Q(:,n-1);Q(:,n-1);Q(:,n-1)];
        Q(1,n) = Q(1,n-1)+dt*b*Qprime(1:2:2*n_a);
        Q(2,n) = Q(2,n-1)+dt*b*Qprime(2:2:2*n_a);
    end
    
    q_ERK4 = Q(1,:);

    plot(t,q_ERK4,'.-g')

    legend_strings = [legend_strings;'ERK4   '];
end

if methods(1) == 7
    % IRK2

    methods(1) = [];
    if isempty(methods)
        methods = 99;
    end

    Q   = zeros(2,N);
    
    Q(:,1) = [1;0];
    
    gamma = 1/2+1/6*sqrt(3);

    a11 = gamma;
    a21 = 1-2*gamma;
    a22 = gamma;
    b1  = 1/2;
    b2  = 1/2;

    A = [ a11  0 
          a21 a22];

    b = [b1 b2];

    n_a = length(A);
    
    AA = zeros(2*n_a);
    for i=1:2*n_a
        if mod(i,2)~=0
            AA(i,1:2:2*n_a) = A(ceil(i/2),:);
        else
            AA(i,2:2:2*n_a) = A(ceil(i/2),:);
        end
    end
    
    RR = zeros(2*n_a);
    for i=1:n_a
        RR(2*i-1:2*i,2*i-1:2*i) = R;
    end
    
    S = inv(inv(RR)-dt*AA);
    
    for n=2:N
        Qprime = S*[Q(:,n-1);Q(:,n-1)];
        Q(1,n) = Q(1,n-1)+dt*b*Qprime(1:2:2*n_a);
        Q(2,n) = Q(2,n-1)+dt*b*Qprime(2:2:2*n_a);
    end

    q_IRK2 = Q(1,:);

    plot(t,q_IRK2,'.-')

    legend_strings = [legend_strings;'IRK2   '];
end

if methods(1) == 8
    % ESDIRK3

    methods(1) = [];
    if isempty(methods)
        methods = 99;
    end

    Q   = zeros(2,N);
    
    Q(:,1) = [1;0];
    
    esdirk3;
    
    n_a = length(rkaI);
    
    AA = zeros(2*n_a);
    for i=1:2*n_a
        if mod(i,2)~=0
            AA(i,1:2:2*n_a) = rkaI(ceil(i/2),:);
        else
            AA(i,2:2:2*n_a) = rkaI(ceil(i/2),:);
        end
    end
    
    RR = zeros(2*n_a);
    for i=1:n_a
        RR(2*i-1:2*i,2*i-1:2*i) = R;
    end
      
    S = inv(eye(2*n_a)-dt*RR*AA);
    
    for n=1:N-1
        Qprime   = S*[Q(:,n);Q(:,n);Q(:,n);Q(:,n)];
        Q(:,n+1) = Qprime(end-1:end);
    end

    q_ESDIRK3 = Q(1,:);

    plot(t,q_ESDIRK3,'.-m')
    
    legend_strings = [legend_strings;'ESDIRK3'];
end

if methods(1) == 9
    % ESDIRK4

    methods(1) = [];
    if isempty(methods)
        methods = 99;
    end

    Q   = zeros(2,N);
    
    Q(:,1) = [1;0];
    
    esdirk4;
    
    n_a = length(rkaI);
    
    AA = zeros(2*n_a);
    for i=1:2*n_a
        if mod(i,2)~=0
            AA(i,1:2:2*n_a) = rkaI(ceil(i/2),:);
        else
            AA(i,2:2:2*n_a) = rkaI(ceil(i/2),:);
        end
    end
    
    RR = zeros(2*n_a);
    for i=1:n_a
        RR(2*i-1:2*i,2*i-1:2*i) = R;
    end
      
    S = inv(eye(2*n_a)-dt*RR*AA);
    
    for n=1:N-1
        Qprime   = S*[Q(:,n);Q(:,n);Q(:,n);Q(:,n);Q(:,n);Q(:,n)];
        Q(:,n+1) = Qprime(end-1:end);
    end

    q_ESDIRK4 = Q(1,:);

    plot(t,q_ESDIRK4,'.-r')

    legend_strings = [legend_strings;'ESDIRK4'];
end

if methods(1) == 10
    % ESDIRK5

    methods(1) = [];
    if isempty(methods)
        methods = 99;
    end

    Q   = zeros(2,N);
    
    Q(:,1) = [1;0];
    
    esdirk5;
    
    n_a = length(rkaI);
    
    AA = zeros(2*n_a);
    for i=1:2*n_a
        if mod(i,2)~=0
            AA(i,1:2:2*n_a) = rkaI(ceil(i/2),:);
        else
            AA(i,2:2:2*n_a) = rkaI(ceil(i/2),:);
        end
    end
    
    RR = zeros(2*n_a);
    for i=1:n_a
        RR(2*i-1:2*i,2*i-1:2*i) = R;
    end
      
    S = inv(eye(2*n_a)-dt*RR*AA);
    
    for n=1:N-1
        Qprime   = S*[Q(:,n);Q(:,n);Q(:,n);Q(:,n);Q(:,n);Q(:,n);Q(:,n);Q(:,n)];
        Q(:,n+1) = Qprime(end-1:end);
    end

    q_ESDIRK5 = Q(1,:);

    plot(t,q_ESDIRK5,'.-c')

    legend_strings = [legend_strings;'ESDIRK5'];
end


legend(legend_strings)