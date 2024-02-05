x1 = linspace(-1,0,50)'*ones(1,100);
y1 = ones(50,1)*linspace(-1,1,100);
y2 = ones(50,1)*linspace(-1,1,100);
x2 = linspace(0,1,50)'*ones(1,100);
x = [x1 ; x2];
y = ones(100,1)*linspace(-1,1,100);
u1 = 1/16*(y1.^2-1).*(x1+1).^2;
u2 = 1/16*(1-y2.^2).*((x2.^2-1)+2*1e-3*x2.*(x2-1));
u = [u1 ; u2];
surf(x,y,u); shading interp