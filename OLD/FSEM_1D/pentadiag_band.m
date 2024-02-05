% Thomas
function x = pentadiag_band(a,b,c,d,e,r)

if length(a)~=length(b)&& length(b)~=length(c) ...
   && length(c)~=length(d) && length(d)~=length(e) && length(e)~=length(r)
    warning('size of vectors a,b,c,d,e & r are not the same'); %#ok<WNTAG>
    return
end

n = length(a);

f = zeros(n,1);
g = zeros(n,1);
s = zeros(n,1);
x = zeros(n,1);

% Downsweep

% first row
f(2) = d(2)/c(1);
g(3) = e(3)/c(1);
s(1) = r(1)/c(1);

den = c(2)-b(1)*f(2);
f(3) = (d(3)-b(1)*g(3))/den;
g(4) = e(4)/den;
s(2) = (r(2)-b(1)*s(1))/den;

for i=3:n-2
   den    = c(i)-a(i-2)*g(i)-(b(i-1)-a(i-2)*f(i-1))*f(i);
   f(i+1) = (d(i+1)-(b(i-1)-a(i-2)*f(i-1))*g(i+1))/den;
   g(i+2) = e(i+2)/den;
   s(i)   = (r(i)-a(i-2)*s(i-2)-(b(i-1)-a(i-2)*f(i-1))*s(i-1))/den;
end
den    = c(n-1)-a(n-3)*g(n-1)-(b(n-2)-a(n-3)*f(n-2))*f(n-1);
f(n)   = (d(n)-(b(n-2)-a(n-3)*f(n-2))*g(n))/den;
s(n-1) = (r(n-1)-a(n-3)*s(n-3)-(b(n-2)-a(n-3)*f(n-2))*s(n-2))/den;
s(n)   = (r(n)-a(n-2)*s(n-2)-(b(n-1)-a(n-2)*f(n-1))*s(n-1))/(c(n)-a(n-2)*g(n)-(b(n-1)-a(n-2)*f(n-1))*f(n));

% Upsweep
x(n) = s(n);
x(n-1) = s(n-1)-f(n)*x(n);
for i=n-2:-1:1
   x(i) = s(i)-f(i+1)*x(i+1)-g(i+2)*x(i+2);
end