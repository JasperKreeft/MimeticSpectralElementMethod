clear all
close all
clc

global N m cc numRows numColumns numElements
global e e_w

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% settings & loop for h-convergence                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e = 0; er = 0;

method = 2;

HconvRange = 2;%[ 2 4 8 16 32 64 ];

NrCellRange = 2;%2:2:14;

m = 2;

% grid option
cc = 0.0;

errorL2_q        = zeros(1,max(length(NrCellRange),length(HconvRange)));
errorL2_q_interp = zeros(1,max(length(NrCellRange),length(HconvRange)));

for Hconv = HconvRange

numRows    = Hconv;
numColumns = Hconv;
numElements = numRows*numColumns;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% code for convergence plots                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(Hconv)

% errorL2        = zeros(1,max(NrCellRange));
% errorL2_interp = zeros(1,max(NrCellRange));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop for p-convergence                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for N = NrCellRange

disp(['N = ' num2str(N)]);

N2 = N*N;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% numbering unknowns                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

localnr = reshape(1:N2,N,N);

globalnr = zeros(numColumns*N,numRows*N);
for j=1:numRows
    for i=1:numColumns
        globalnr((i-1)*N+(1:N),(j-1)*N+(1:N)) = ((i-1)+(j-1)*numColumns)*N2+localnr;
    end
end

% disp('  '); disp('globalnr = '); disp('  '); disp(num2str(flipud(globalnr')));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mesh generation                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xiGLL,wGLL] = GLLnodes(N);    etaGLL = xiGLL;  % Gauss-Lobotto-Legendre
[xiG,wG]   = Gnodes(N);        etaG   = xiG;    % Gauss
xiEG     = [-1 xiG 1];         etaEG  = xiEG;   % Extended Gauss

XiGLLGLL = xiGLL'*ones(1,N+1); EtaGLLGLL = XiGLLGLL';
XiGLLG  = xiGLL'*ones(1,N);    EtaGLLG   = ones(N+1,1)*etaG;
XiGGLL  = xiG'*ones(1,N+1);    EtaGGLL   = ones(N,1)*etaGLL;
XiEGEG  = xiEG'*ones(1,N+2);   EtaEGEG   = XiEGEG';

[h  ,e  ] = MimeticpolyVal(xiGLL,N,1);
[h_w,e_w] = MimeticpolyVal(xiEG ,N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('assembly matrices and rhs')

% Main mesh creation. Here for multi-element on unit square [0,1]x[0,1] or
% standard square [-1,1]x[-1,1]

xibLR  = linspace(-1,1,numColumns+1);   % This might be more advanced !!!!!
etabAB = linspace(-1,1,numRows+1)   ;   % This might be more advanced !!!!!

% Global mesh
XGLLGLL = zeros(N*numColumns+1,N*numRows+1);
YGLLGLL = zeros(N*numColumns+1,N*numRows+1);
XGG   = zeros(N*numColumns  ,N*numRows  );
YGG   = zeros(N*numColumns  ,N*numRows  );
ind1 = 2*N*(N+1)*numColumns*numRows;
ind2 = 4*N^2*(N+1)^2*numColumns*numRows;
if method == 1
    H = spalloc(ind1,ind1,ind2);
elseif method == 2
    B = spalloc(ind1,ind1,ind2);
    ind1 = (N2+3*N)*numColumns*numRows;
    ind2 = (N2+3*N)^2*numColumns*numRows;
    C = spalloc(ind1,ind1,ind2);
end
ind = numColumns*numRows*(N2+2*N)+(numRows+numColumns)*N;
f = zeros(ind,1);
for r=1:numRows
    for c=1:numColumns

        XibGLLGLL  = ( (xibLR(c)+xibLR(c+1))/2+(xibLR(c+1)-xibLR(c))/2*xiGLL )'*ones(1,N+1);
        EtabGLLGLL = ones(N+1,1)*( (etabAB(r)+etabAB(r+1))/2+(etabAB(r+1)-etabAB(r))/2*etaGLL );

        XibGG  = ( (xibLR(c)+xibLR(c+1))/2+(xibLR(c+1)-xibLR(c))/2*xiG )'*ones(1,N);
        EtabGG = ones(N,1)*( (etabAB(r)+etabAB(r+1))/2+(etabAB(r+1)-etabAB(r))/2*etaG );

        % gridtype = 'sinecurve' % This might be more advanced !!!!!
        XeGLLGLL = XibGLLGLL + cc*sin(pi*XibGLLGLL).*sin(pi*EtabGLLGLL);
        YeGLLGLL = EtabGLLGLL+ cc*sin(pi*XibGLLGLL).*sin(pi*EtabGLLGLL);

        XeGG = XibGG + cc*sin(pi*XibGG).*sin(pi*EtabGG);
        YeGG = EtabGG+ cc*sin(pi*XibGG).*sin(pi*EtabGG);

        if c<numColumns && r<numRows
            XGLLGLL((c-1)*N+(1:N),(r-1)*N+(1:N)) = XeGLLGLL(1:N,1:N);
            YGLLGLL((c-1)*N+(1:N),(r-1)*N+(1:N)) = YeGLLGLL(1:N,1:N);
        elseif c<numColumns && r==numRows
            XGLLGLL((c-1)*N+(1:N),(r-1)*N+(1:N+1)) = XeGLLGLL(1:N,1:N+1);
            YGLLGLL((c-1)*N+(1:N),(r-1)*N+(1:N+1)) = YeGLLGLL(1:N,1:N+1);
        elseif c==numColumns && r<numRows
            XGLLGLL((c-1)*N+(1:N+1),(r-1)*N+(1:N)) = XeGLLGLL(1:N+1,1:N);
            YGLLGLL((c-1)*N+(1:N+1),(r-1)*N+(1:N)) = YeGLLGLL(1:N+1,1:N);
        elseif c==numColumns && r==numRows
            XGLLGLL((c-1)*N+(1:N+1),(r-1)*N+(1:N+1)) = XeGLLGLL;
            YGLLGLL((c-1)*N+(1:N+1),(r-1)*N+(1:N+1)) = YeGLLGLL;
        end
        XGG((c-1)*N+(1:N),(r-1)*N+(1:N)) = XeGG;
        YGG((c-1)*N+(1:N),(r-1)*N+(1:N)) = YeGG;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Metric relations                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if method == 1
            He = elementmatrix_method1(N,xiGLL,xiEG,wG,wGLL,XibGLLGLL,EtabGLLGLL);
            ind = 2*N*(N+1)*((r-1)*numColumns+c-1)+(1:2*N*(N+1));
            H(ind,ind) = He;

        elseif method == 2
            [Be,Ce] = elementmatrix_method2(N,xiGLL,xiEG,wG,wGLL,XibGLLGLL,EtabGLLGLL);
            ind = 2*N*(N+1)*((r-1)*numColumns+c-1)+(1:2*N*(N+1));
            B(ind,ind) = Be;
            
            ind = (c-1+(r-1)*numColumns)*(N2+4*N)+(1:N2+4*N);
            C(ind,ind) = Ce;

        end


        % keyboard
        % Method 2 werkt nog niet, omdat He niet bestaat!!


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forcing function                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        F = force(XibGLLGLL(:,1),EtabGLLGLL(1,:));
        i = (c-1)*N+(1:N);
        j = (r-1)*N+(1:N);
        ind = globalnr(i,j)+...
              +(r-1)*(2*numColumns+1)*N+...
              +(r==1)*(c>1)*(c-1)*N+(r>1)*numColumns*N+...
              +(c-1)*2*N+(c>1)*N;
        f(ind) = F;

    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nodes_in_element = N2+4*N;
edges_in_element = 2*N*(N+1);
cells_in_element = N2;

if method == 2

for r = numRows:-1:1
    for c = numColumns:-1:1
        max_in_element       = (c+(r-1)*numColumns)*nodes_in_element;
        % max_in_left_element  = ((c-1)+(r-1)*numColumns)*nodes_in_element;
        % max_in_lower_element = (c+(r-2)*numColumns)*nodes_in_element;

        if r>1
%             for i=1:N
                ind2 = max_in_element-(2*(1:N)-1);
                C(:,ind2) = [];
                C(ind2,:) = [];
%             end
        end
        if c>1
%             for i=1:N
                ind2 = max_in_element-2*N-(2*(1:N)-1);
                C(:,ind2) = [];
                C(ind2,:) = [];
%             end
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Topology relations                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('assembly topology')

% filename = ['Dp_H' num2str(Hconv) '_N' num2str(N) '.mat'];
% load(filename);


Dpe = topology(N);

% Divergence operator

Dp = kron(speye(numRows*numColumns),Dpe);

for r = numRows:-1:1
    for c = numColumns:-1:1
        max_in_element       = (c+(r-1)*numColumns)*nodes_in_element;
        max_in_left_element  = ((c-1)+(r-1)*numColumns)*nodes_in_element;
        max_in_lower_element = (c+(r-2)*numColumns)*nodes_in_element;

        if r>1
            for i=1:N
                ind1 = max_in_lower_element-2*(i-1);
                ind2 = max_in_element-(2*i-1);
                Dp(ind1,:) = Dp(ind1,:) + Dp(ind2,:);
                Dp(ind2,:) = [];
            end            
        end
        if c>1
            for i=1:N
                ind1 = max_in_left_element-2*N-2*(i-1);
                ind2 = max_in_element-2*N-(2*i-1);
                Dp(ind1,:) = Dp(ind1,:) + Dp(ind2,:);
                Dp(ind2,:) = [];
            end
        end
    end
end

% Gradient operator
Gd = -Dp';
% 
% save(['Dp_H' num2str(Hconv) '_N' num2str(N)],'Dp');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct matrix                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('solve part 1')
pause(0.05)
if method == 1
    Q = H*Gd;        clear H;
    A = Dp*Q;          clear Q;
elseif method == 2
    A = -Dp*(B\Dp');% clear B;
%     A = -Dp*inv(B)*Dp'; clear B;
end
keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boundary conditions                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Hconv>1 %&& method == 1
% remove bc columns

max_in_element = (numRows*numColumns)*N2+(numRows*(numColumns-1)+numColumns*(numRows-1))*N+... % internal points
                 2*N*(numColumns+numRows); % boundary points
for r = numRows:-1:1
    for c = numColumns:-1:1
        if r == numRows
            for i=N:-1:1
                ind = max_in_element-(N-i)*(1+(r==1));
                A( : , ind ) = [];
                A( ind , : ) = [];
%                 C( : , ind ) = [];
%                 C( ind , : ) = [];
                f(ind) = [];
            end
        end
        if r == 1
            for i=N:-1:1
                ind = max_in_element-(1+2*(N-i))-(r==numRows)*(i-1);
                A( : , ind ) = [];
                A( ind , : ) = [];
%                 C( : , ind ) = [];
%                 C( ind , : ) = [];
                f(ind) = [];
            end
        end
        if c == numColumns
            for i=N:-1:1
                ind = max_in_element-N-(N-i)-N*(r==1)-(N-i)*(c==1);
                A( : , ind ) = [];
                A( ind , : ) = [];
%                 C( : , ind ) = [];
%                 C( ind , : ) = [];
                f(ind) = [];
            end
        end
        if c == 1
            for i=N:-1:1
                ind = max_in_element-N-1-2*(N-i)-N*(r==1)-(i-1)*(c==numColumns);
                A( : , ind ) = [];
                A( ind , : ) = [];
%                 C( : , ind ) = [];
%                 C( ind , : ) = [];
                f(ind) = [];
            end
        end
        max_in_element = max_in_element-(N2+2*N)-N*((c==1)+(r==1));
    end
end
end

%%%%%%%%%%%%%%%%%%%%%

% fbc = zeros((numRows*(numColumns-1)+numColumns*(numRows-1))*N,1);
% 
% f = [ f ; fbc];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solving the system                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi = A\f;

if method == 2
    
    phi_in = phi;
    phi = zeros(numRows*numColumns*(N2+2*N)+(numRows+numColumns)*N,1);
    
    inner = [];
    for r=1:numRows
        for c=1:numColumns
            max_prev_element = ((r-1)*numColumns+c-1)*(N2+2*N)+((r==1)*(c-1)+(r>1)*numColumns)*N;
            inner = [inner max_prev_element+(r-1+(c>1))*N+(1:N2)];
            if c<numColumns
                inner = [inner max_prev_element+(r-1+(c>1))*N+N2+((1+(c==1)):(1+(c==1)):(1+(c==1))*N)];
            end
            if r<numRows
                inner = [inner max_prev_element+(r-1+(c>1))*N+N2+N+(c==1)*N+((1+(r==1)):(1+(r==1)):(1+(r==1))*N)];
            end
        end
    end

    phi(inner) = phi_in;
    
%     Q = -(B\Dp')*phi;
    
    disp('solve part 2')
    
    phi = C\phi;

    phi = phi(inner);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Postprocessen                                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('post-processen')

postprocessen_me

% postprocessen_q

end % for N
end % for H

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% code for convergence plots                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% errorL2

figure(1)

if length(NrCellRange)>1
    semilogy(NrCellRange,errorL2)
    hold on
    semilogy(NrCellRange,errorL2_interp,'--r')
    filename = ['Pconv_H' num2str(Hconv) '_c' num2str(cc) '.mat'];
end
if length(HconvRange)>1
    loglog(2./(HconvRange),errorL2)
    hold on
    loglog(2./(HconvRange),errorL2_interp,'--xr')
    filename = ['Hconv_N' num2str(N) '_c' num2str(cc) '.mat'];
end
errorL2

errorL2_interp

% figure(2)
% 
% if length(NrCellRange)>1
%     semilogy(NrCellRange,errorL2_q)
%     hold on
%     semilogy(NrCellRange,errorL2_q_interp,'--r')
% end
% if length(HconvRange)>1
%     loglog(2./(HconvRange),errorL2_q)
%     hold on
%     loglog(2./(HconvRange),errorL2_q_interp,'--xr')
% end
% errorL2_q
% 
% errorL2_q_interp


save(filename,'N','NrCellRange','errorL2','errorL2_interp','Hconv');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
