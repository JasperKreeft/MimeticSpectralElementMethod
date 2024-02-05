function [qphi_bc,q_bc,ind_bc,sindq]=boundaryconditions_ConvectionVolume(globalnr_1h,globalnr_1v,XGLLGLL,YGLLGLL,V)

global N numRows numColumns

a = V(1);
b = V(2);

% indices

% Boundary indices for matrix A (all k-cells numbered)
ind_leftbc  = globalnr_1h(1,:)   ; ind_leftbc(ind_leftbc==0)   = [];
ind_lowerbc = globalnr_1v(:,1)'  ; ind_lowerbc(ind_lowerbc==0) = [];
ind_q_bc = sort([ind_leftbc ind_lowerbc]);

ind_bc = ind_q_bc;

% Boundary indices for qphi_bc (only k-cells numbered with boundary condition)
ind_bc_left  = [(1:N) reshape((ones(numRows-1,1)*(1:N)+(0:N:N*(numRows-2))'*ones(1,N)+N*(numColumns+1))',1,[])];
ind_bc_lower = [N+1:N*numColumns N*numColumns+(1:N)];
indq   = [ind_bc_left ind_bc_lower]; indq(indq == 0) = [];
sindq  = size(indq,2);

%%%

q_leftbc  = numRows*a/pi*(sin(pi*(1+YGLLGLL(1,2:end)))-sin(pi*(1+YGLLGLL(1,1:end-1))));
q_lowerbc = numColumns*b/pi*(sin(pi*(XGLLGLL(2:end,1)+1))-sin(pi*(XGLLGLL(1:end-1,1)+1)));

% q_leftbc  = zeros(1,numRows*N);
% q_lowerbc = zeros(numColumns*N,1);
% for i=1:numColumns*N
% if XGLLGLL(i,1)>0
%     q_lowerbc(i,1) = 0;
% elseif XGLLGLL(i+1,1)>0
%     q_lowerbc(i,1) = -XGLLGLL(i,1)+1/(2*pi)*sin(2*pi*XGLLGLL(i,1));
% else
%     q_lowerbc(i,1) = (XGLLGLL(i+1,1)-XGLLGLL(i,1))-1/(2*pi)*(sin(2*pi*XGLLGLL(i+1,1))-sin(2*pi*XGLLGLL(i,1)));
% end
% end


% q_leftbc  = zeros(1,numRows*N);
% q_lowerbc = zeros(numColumns*N,1);
% for i=1:numColumns*N
%     xstep = -0.4;
%     % step-function
%     if XGLLGLL(i+1,1)<=xstep
%         q_lowerbc(i,1) = 0;
%     elseif XGLLGLL(i,1)<=xstep && XGLLGLL(i+1,1)>xstep
%         q_lowerbc(i,1) = XGLLGLL(i+1,1)-xstep;
%     elseif XGLLGLL(i,1)>xstep
%         q_lowerbc(i,1) = XGLLGLL(i+1,1)-XGLLGLL(i,1);
%     end
% end

% q_leftbc  = numRows*a*diff(YGLLGLL(1,:));
% q_lowerbc = zeros(numColumns*N,1);

% keyboard
q_bc = [q_leftbc' ; q_lowerbc];
qphi_bc(indq,1) = q_bc;

%%%