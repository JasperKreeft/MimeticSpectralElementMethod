% Thomas
function x = tridiag_band(a,b,c,r)

if length(a)~=length(b)&& length(a)~=length(c) && length(c)~=length(r)
    warning('size of vectors a,b,c & r are not the same'); %#ok<WNTAG>
    return
end

n = length(a);

e = zeros(n,1);
s = zeros(n,1);
x = zeros(n,1);

% Downsweep
e(1) = 1;
s(1) = r(1)/b(1);
e(2) = c(2)/b(1);
for i=2:n-1
   den = b(i)-a(i-1)*e(i);
   s(i) = (r(i)-a(i-1)*s(i-1))/den;
   e(i+1) = c(i+1)/den;
end
s(n) = (r(n)-a(n-1)*s(n-1))/(b(n)-a(n-1)*e(n));


% Upsweep
x(n) = s(n);
for i=n-1:-1:1
   x(i) = s(i)-e(i+1)*x(i+1);
end