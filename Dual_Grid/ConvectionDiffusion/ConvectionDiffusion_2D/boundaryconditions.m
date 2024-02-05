function [qphi_bc,q_bc,phi_bc,ind_bc,sindq]=boundaryconditions(globalnr_0,globalnr_1h,globalnr_1v,XGLLGLL,YGLLGLL,XGG,YGG)

global N numRows numColumns RC
global edges_in_element
global BC

%      le ri lo up
BC = [ 1  1  1  1      % q
       0  0  0  0 ];   % phi

% indices

% Boundary indices for matrix A (all k-cells numbered)
ind_leftbc  = BC(1,1)*globalnr_1h(1,:)   ; ind_leftbc(ind_leftbc==0)   = [];
ind_rightbc = BC(1,2)*globalnr_1h(end,:) ; ind_rightbc(ind_rightbc==0) = [];
ind_lowerbc = BC(1,3)*globalnr_1v(:,1)'  ; ind_lowerbc(ind_lowerbc==0) = [];
ind_upperbc = BC(1,4)*globalnr_1v(:,end)'; ind_upperbc(ind_upperbc==0) = [];
ind_q_bc = sort([ind_leftbc ind_rightbc ind_lowerbc ind_upperbc]);

ind_leftbc  = BC(2,1)*globalnr_0(1,:)   ; ind_leftbc(ind_leftbc==0)   = [];
ind_rightbc = BC(2,2)*globalnr_0(end,:) ; ind_rightbc(ind_rightbc==0) = [];
ind_lowerbc = BC(2,3)*globalnr_0(:,1)'  ; ind_lowerbc(ind_lowerbc==0) = [];
ind_upperbc = BC(2,4)*globalnr_0(:,end)'; ind_upperbc(ind_upperbc==0) = [];
ind_phi_bc = sort([ind_leftbc ind_rightbc ind_lowerbc ind_upperbc])+RC*edges_in_element;

ind_bc = [ind_q_bc ind_phi_bc];

% Boundary indices for qphi_bc (only k-cells numbered with boundary condition)
ind_bc_left  = BC(1,1)*[(1:N) reshape((ones(numRows-1,1)*(1:N)+(0:N*(1+BC(1,2)):(1+BC(1,2))*N*(numRows-2))'*ones(1,N)+N*(BC(1,3)*numColumns+BC(1,2)+1))',1,[])];
ind_bc_right = BC(1,2)*[N*(BC(1,1)+BC(1,3)*(numColumns-1))+(1:N) reshape((ones(numRows-2,1)*(1:N)+((BC(1,1)*N):(BC(1,1)+1)*N:(BC(1,1)+1)*N*(numRows-3+BC(1,1)))'*ones(1,N)+N*(BC(1,1)+BC(1,3)*numColumns+1))',1,[]) N*(BC(1,1)*numRows+(numRows-1)+BC(1,3)*numColumns+BC(1,4)*(numColumns-1))+(1:N)];
ind_bc_lower = BC(1,3)*[BC(1,1)*N+1:N*(numColumns-1+BC(1,1)) N*(numColumns-1+BC(1,1)+BC(1,2))+(1:N)];
ind_bc_upper = BC(1,4)*[N*(BC(1,1)*numRows+BC(1,2)*(numRows-1)+BC(1,3)*numColumns)+(1:N*(numColumns-1)) N*((BC(1,1)+BC(1,2))*numRows+BC(1,3)*numColumns+numColumns-1)+(1:N)];

indq   = [ind_bc_left ind_bc_right ind_bc_lower ind_bc_upper]; indq(indq == 0)   = [];
sindq  = size(indq,2);

ind_bc_left  = BC(2,1)*[(1:N) reshape((ones(numRows-1,1)*(1:N)+(0:N*(1+BC(2,2)):(1+BC(2,2))*N*(numRows-2))'*ones(1,N)+N*(BC(2,3)*numColumns+BC(2,2)+1))',1,[])];
ind_bc_right = BC(2,2)*[N*(BC(2,1)+BC(2,3)*(numColumns-1))+(1:N) reshape((ones(numRows-2,1)*(1:N)+((BC(2,1)*N):(BC(2,1)+1)*N:(BC(2,1)+1)*N*(numRows-3+BC(2,1)))'*ones(1,N)+N*(BC(2,1)+BC(2,3)*numColumns+1))',1,[]) N*(BC(2,1)*numRows+(numRows-1)+BC(2,3)*numColumns+BC(2,4)*(numColumns-1))+(1:N)];
ind_bc_lower = BC(2,3)*[BC(2,1)*N+1:N*(numColumns-1+BC(2,1)) N*(numColumns-1+BC(2,1)+BC(2,2))+(1:N)];
ind_bc_upper = BC(2,4)*[N*(BC(2,1)*numRows+BC(2,2)*(numRows-1)+BC(2,3)*numColumns)+(1:N*(numColumns-1)) N*((BC(2,1)+BC(2,2))*numRows+BC(2,3)*numColumns+numColumns-1)+(1:N)];

indphi = [ind_bc_left ind_bc_right ind_bc_lower ind_bc_upper]; indphi(indphi==0) = [];
sindphi = size(indphi,2); %#ok<NASGU>

%%%

if BC(1,1)==true
    q_leftbc = (cos(pi*XGLLGLL(1,1:end-1)).*(cos(pi*YGLLGLL(1,1:end-1))-cos(pi*YGLLGLL(1,2:end))))';
else
    q_leftbc = [];
end
if BC(1,2)==true
    q_rightbc = (cos(pi*XGLLGLL(end,1:end-1)).*(cos(pi*YGLLGLL(end,1:end-1))-cos(pi*YGLLGLL(end,2:end))))';
else
    q_rightbc = [];
end
if BC(1,3)==true
    q_lowerbc = (cos(pi*XGLLGLL(1:end-1,1))-cos(pi*XGLLGLL(2:end,1))).*cos(pi*YGLLGLL(1:end-1,1));%        0*YGLLGLL(1:end-1,1);
else
    q_lowerbc = [];
end
if BC(1,4)==true
    q_upperbc = (cos(pi*XGLLGLL(1:end-1,end))-cos(pi*XGLLGLL(2:end,end))).*cos(pi*YGLLGLL(1:end-1,end)); %0*YGLLGLL(1:end-1,end);
else
    q_upperbc = [];
end
q_bc = [q_leftbc ; q_rightbc ; q_lowerbc ; q_upperbc];
qphi_bc(indq,1) = q_bc;

%%%

if BC(2,1)==true
    phi_leftbc =0*(1+YGG(1,:)')/2;% 0.5*ones(size(YGG(1,:)')); %
else
    phi_leftbc = [];
end
if BC(2,2)==true
    phi_rightbc = 0*(1+YGG(end,:)')/2;%0.5*ones(size(YGG(end,:)')); %
else
    phi_rightbc = [];
end
if BC(2,3)==true
    phi_lowerbc = 0*0*XGG(:,1); % 1-XGG(:,1).^2;%0.5*ones(size(XGG(:,1)));%
else
    phi_lowerbc = [];
end
if BC(2,4)==true
    phi_upperbc = 0*ones(size(XGG(:,end)));%*XGG(:,end);
else
    phi_upperbc = [];
end
phi_bc = [phi_leftbc ; phi_rightbc ; phi_lowerbc ; phi_upperbc];
qphi_bc(indphi+sindq,1) = phi_bc;

%%%%