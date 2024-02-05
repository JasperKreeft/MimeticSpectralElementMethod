function [k11 k12 k21 k22 dk11 dk12 dk21 dk22] = Kmatrix(X,Y)


% k11 = 10*ones(size(X));
% k12 = 3*ones(size(X));
% k21 = 3*ones(size(X));
% k22 = ones(size(X));
% 
% dk11 = zeros(size(X));
% dk12 = zeros(size(X));
% dk21 = zeros(size(X));
% dk22 = zeros(size(X));


% k11 = X+eps;
% k12 = 0*ones(size(X));
% k21 = 0*ones(size(X));
% k22 = ones(size(X));
% 
% dk11 = ones(size(X));
% dk12 = zeros(size(X));
% dk21 = zeros(size(X));
% dk22 = zeros(size(X));


k11 = (X+1).^2+Y.^2;
k12 = sin(X.*Y);
k21 = sin(X.*Y);
k22 = (X+1).^2;

dk11 = [];
dk12 = [];
dk21 = [];
dk22 = [];



% global i
% 
% part1 = []; part2 = [];
% for j=1:3
% 
% part1 = [ part1 [1:3  7:9  13:15]+(j-1)*15 ];
% part1 = [ part1     [ 4:6 10:12 ]+(j+2)*15 ];
% part1 = [ part1 [1:3  7:9  13:15]+(j+5)*15 ];
% part1 = [ part1     [ 4:6 10:12 ]+(j+8)*15 ];
% part1 = [ part1 [1:3  7:9  13:15]+(j+11)*15 ];
% part2 = [ part2     [ 4:6 10:12 ]+(j-1)*15 ];
% part2 = [ part2 [1:3  7:9  13:15]+(j+2)*15 ];
% part2 = [ part2     [ 4:6 10:12 ]+(j+5)*15 ];
% part2 = [ part2 [1:3  7:9  13:15]+(j+8)*15 ];
% part2 = [ part2     [ 4:6 10:12 ]+(j+11)*15 ];
% 
% end
% 
% part1 = sort(part1);
% part2 = sort(part2);
% 
% k12 = zeros(size(X)); k21 = k12;
% 
% if sum(part1==i)
% k11 = 1*ones(size(X));
% elseif sum(part2==i)
% k11 = 2;
% end
% k22 = k11;
% 
% dk11 = [];
% dk12 = [];
% dk21 = [];
% dk22 = [];
% 
