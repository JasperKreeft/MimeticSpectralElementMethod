function [qphi_bc,ind_phi_bc2]=boundaryconditions(globalnr_0,XGLLGLL,YGLLGLL)

global N numRows numColumns

ind_leftbc  = globalnr_0(1,:)   ;
ind_rightbc = globalnr_0(end,:) ;
ind_lowerbc = globalnr_0(:,1)'  ;
ind_upperbc = globalnr_0(:,end)';
ind_phi_bc2 = sort([ ind_leftbc ind_rightbc ind_lowerbc ind_upperbc]); ind_phi_bc2(ind_phi_bc2==0) = [];

% Boundary indices for qphi_bc (only k-cells numbered with boundary condition)
ind_bc_left  = [(1:N) reshape((ones(numRows-1,1)*(1:N)+(0:2*N:2*N*(numRows-2))'*ones(1,N)+N*(numColumns+2))',1,[])];
ind_bc_right = [N*(1+(numColumns-1))+(1:N) reshape((ones(numRows-2,1)*(1:N)+(N:2*N:2*N*(numRows-2))'*ones(1,N)+N*(numColumns+2))',1,[]) N*(numRows+(numRows-1)+numColumns+(numColumns-1))+(1:N)];
ind_bc_lower = [N+1:N*numColumns N*(numColumns+1)+(1:N)];
ind_bc_upper = [N*(numRows+(numRows-1)+numColumns)+(1:N*(numColumns-1)) N*(2*numRows+numColumns+numColumns-1)+(1:N)];

indq   = [ind_bc_left ind_bc_right ind_bc_lower ind_bc_upper];

%%%

q_leftbc = (cos(pi*XGLLGLL(1,1:end-1)).*(cos(pi*YGLLGLL(1,1:end-1))-cos(pi*YGLLGLL(1,2:end))))';
q_rightbc = (cos(pi*XGLLGLL(end,1:end-1)).*(cos(pi*YGLLGLL(end,1:end-1))-cos(pi*YGLLGLL(end,2:end))))';
q_lowerbc = (cos(pi*XGLLGLL(1:end-1,1))-cos(pi*XGLLGLL(2:end,1))).*cos(pi*YGLLGLL(1:end-1,1));
q_upperbc = (cos(pi*XGLLGLL(1:end-1,end))-cos(pi*XGLLGLL(2:end,end))).*cos(pi*YGLLGLL(1:end-1,end));

q_bc = [q_leftbc ; -q_rightbc ; q_lowerbc ; -q_upperbc];
qphi_bc(indq,1) = q_bc;

%%%