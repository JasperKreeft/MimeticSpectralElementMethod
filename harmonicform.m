clear all
close all
clc

E10 = [ -1,1,0,0,0,0,0,0
        -1,0,0,0,0,0,1,0
        -1,0,1,0,0,0,0,0
        0,-1,0,1,0,0,0,0
        0,-1,0,0,0,0,0,1
        0,0,-1,1,0,0,0,0
        0,0,-1,0,1,0,0,0
        0,0,0,-1,0,1,0,0
        0,0,0,0,-1,1,0,0
        0,0,0,0,-1,0,1,0
        0,0,0,0,0,-1,0,1
        0,0,0,0,0,0,-1,1 ];

E21 = [ 1, 0,-1, 1,0,-1,0, 0,0, 0, 0, 0
        0,-1, 1, 0,0, 0,1, 0,0, 1, 0, 0
        0, 0, 0,-1,1, 0,0,-1,0, 0,-1, 0
        0, 0, 0, 0,0, 0,0, 0,1,-1, 1,-1 ];
    
E01 = E10';
E12 = E21';

f = null(E21);


% E = E10'*E10;
% for n=1:8
%     for k=1:8 
%         fe(k,1) = dot(f(:,n),E10(:,k));
%     end
%     a(:,n) = E\fe;
%     
%     h(:,n) = f(:,n)-
%     
% end

n=1;
k=1;
fe = dot(f(:,n),E10(:,k));
for i=1:8
    ee(i,1) = dot(E10(:,i),E10(:,k));
end
a = ee/fe;


h=f(:,n);
for i=1:8
    h=h-a(i)*E10(:,i);
end
h


n=2;
k=1;
fe = dot(f(:,n),E10(:,k));
for i=1:8
    ee(i,1) = dot(E10(:,i),E10(:,k));
end
a = ee/fe;


h=f(:,n);
for i=1:8
    h=h-a(i)*E10(:,i);
end
h