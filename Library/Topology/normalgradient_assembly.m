function NG = normalgradient_assembly()

global N nr_0 nr_1
global numElements
global globalnr_0 globalnr_1v globalnr_1h

disp('assembly normal gradient')

NGe = normalgrad(N);

[sr,sc,nge]=find(NGe);

ng = zeros(1,2*numElements*2*N*(N+1));
spind01_c = ng;
spind01_r = ng;

for i=1:numElements

    ind_01 = (1:4*N*(N+1)) + (i-1)*4*N*(N+1);

    spind01_c(ind_01) = globalnr_0(sc,i);
    globalnr_1 = [ globalnr_1v ; globalnr_1h ];
    spind01_r(ind_01) = globalnr_1(sr,i);
    ng(ind_01) = nge;

end

NG = sparse(spind01_r,spind01_c,ng,nr_1,nr_0);

NG = sign(NG);

%% full-memory version
% NG = zeros(nr_1,nr_0);
% 
% for i=1:numElements
% 
%     ind0 = globalnr_0(:,i);
%     ind1 = [ globalnr_1v(:,i) ; globalnr_1h(:,i) ];
% 
%     % Gradient operator
%     NG(ind1,ind0) = NGe;
% 
% end
% 
% NG = sparse(NG);


%% Trage versie
% 
% [sr,sc,nge]=find(NGe);
% 
% ng = zeros(1,2*numElements*2*N*(N+1));
% spind01_c = ng;
% spind01_r = ng;
% 
% for i=1:numElements
% 
%     ind_01 = (1:4*N*(N+1)) + (i-1)*4*N*(N+1);
% 
%     spind01_c(ind_01) = globalnr_0(sc,i);
%     globalnr_1 = [ globalnr_1v ; globalnr_1h ];
%     spind01_r(ind_01) = globalnr_1(sr,i);
%     ng(ind_01) = nge;
% 
% end
% disp('check')
% for k=1:nr_0
%     for l=1:nr_1
%         test = spind01_c==k & spind01_r==l;
%         if sum(test)>1
%             [~,c] = find(test);
%             spind01_c(c(2:end)) = [];
%             spind01_r(c(2:end)) = [];
%             ng(c(2:end))        = [];
%         end
%     end
% end
% disp('checkcheck')
% NG = sparse(spind01_r,spind01_c,ng,nr_1,nr_0);
