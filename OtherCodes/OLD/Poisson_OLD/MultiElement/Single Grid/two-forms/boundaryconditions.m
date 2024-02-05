function [pu_bc,ind_u_bc,ind_p_bc] = boundaryconditions(globalnr_1h,globalnr_1v,globalnr_2,X,Y)

global N m numRows numColumns

ind_leftbc  = globalnr_1v(1,:)   ;
ind_rightbc = globalnr_1v(end,:) ;
ind_lowerbc = globalnr_1h(:,1)'  ;
ind_upperbc = globalnr_1h(:,end)';
ind_p_bc = [ ind_leftbc ind_rightbc ind_lowerbc ind_upperbc ];


ind_leftbc  = globalnr_2(1,:)   ;
ind_rightbc = globalnr_2(end,:) ;
ind_lowerbc = globalnr_2(:,1)'  ;
ind_upperbc = globalnr_2(:,end)';
ind_u_bc = [ ind_leftbc ind_rightbc ind_lowerbc ind_upperbc];

% Boundary indices for qphi_bc (only k-cells numbered with boundary condition)
ind_bc_left  = [(1:N) reshape((ones(numRows-1,1)*(1:N)+(0:2*N:2*N*(numRows-2))'*ones(1,N)+N*(numColumns+2))',1,[])];
ind_bc_right = [N*(1+(numColumns-1))+(1:N) reshape((ones(numRows-2,1)*(1:N)+(N:2*N:2*N*(numRows-2))'*ones(1,N)+N*(numColumns+2))',1,[]) N*(numRows+(numRows-1)+numColumns+(numColumns-1))+(1:N)];
ind_bc_lower = [N+1:N*numColumns N*(numColumns+1)+(1:N)];
ind_bc_upper = [N*(numRows+(numRows-1)+numColumns)+(1:N*(numColumns-1)) N*(2*numRows+numColumns+numColumns-1)+(1:N)];

indp   = [ind_bc_left ind_bc_right ind_bc_lower ind_bc_upper];

p_leftbc = (-1)^(m+1)*(cos(m*pi*Y(1,2:end))-cos(m*pi*Y(end,1:end-1)))';
p_rightbc = (-1)^(m+1)*(cos(m*pi*Y(1,2:end))-cos(m*pi*Y(end,1:end-1)))';
p_lowerbc = (-1)^(m+1)*(cos(m*pi*X(2:end,1))-cos(m*pi*X(1:end-1,end)));
p_upperbc = (-1)^(m+1)*(cos(m*pi*X(2:end,1))-cos(m*pi*X(1:end-1,end)));

p_bc = [p_leftbc ; -p_rightbc ; p_lowerbc ; -p_upperbc];
pu_bc(indp,1) = p_bc;

%%%