function BC = boundary_oneforms(r,c,FunctionType,Domain,DomInfo,boundary,OutputSelect)

global N numElements numRows numColumns
global xi
global n xibar Gw

if ~exist('numElements','var')
    numElements = 1;
end
if isempty(numRows)
numRows = sqrt(numElements);
numColumns = numRows;
end

eta = xi;
etabar = xibar;

% xglobal = linspace(-1,1,numcolumns+1);
% yglobal = linspace(-1,1,numRows+1)   ;

% Main Mesh creation. Here for multi-element on unit square [0,1]x[0,1] or
% standard square [-1,1]x[-1,1]

xibLR  = kron(ones(1,numRows+1),linspace(-1,1,numColumns+1));
etabAB = kron(linspace(-1,1,numRows+1),ones(1,numColumns+1));

ind = c+(r-1)*(numColumns+1);
xiLR  = (xibLR(ind)+xibLR(ind+1))/2+(xibLR(ind+1)-xibLR(ind))/2*xi';
etaAB = (etabAB(ind)+etabAB(ind+numColumns+1))/2+(etabAB(ind+numColumns+1)-etabAB(ind))/2*eta;


arg1 = strcmp(boundary,'left');
arg2 = strcmp(boundary,'right');
arg3 = strcmp(boundary,'below');
arg4 = strcmp(boundary,'above');

switch OutputSelect

    case {'zero','two'} 

        if arg1 || arg2
            Uxi = zeros(N,1);
            for j = 1:N
                Xi  =  ((-1)*arg1+(+1)*arg2)*ones(n,1);
                Eta = (etaAB(j)+etaAB(j+1))/2+(etaAB(j+1)-etaAB(j))/2*etabar;

                [X,Y] = coordinatemap(Xi,Eta,Domain,DomInfo);

                u = exact_solution(X,Y,FunctionType,OutputSelect);
                u = u + exact_solution(X,Y,FunctionType,'mass'); disp('Warning: Mass added')

                Jbar = (eta(j+1)-eta(j))/2;
                Uxi(j,1) = sum( u.*Gw*Jbar );
                BC = Uxi;
            end

        elseif arg3 || arg4
            Ueta = zeros(N,1);
            for i = 1:N
                Xi  = (xiLR(i)+xiLR(i+1))/2+(xiLR(i+1)-xiLR(i))/2*xibar;
                Eta =  ((-1)*arg3+(+1)*arg4)*ones(n,1);

                [X,Y] = coordinatemap(Xi,Eta,Domain,DomInfo);

                u = exact_solution(X,Y,FunctionType,OutputSelect);
                u = u + exact_solution(X,Y,FunctionType,'mass'); disp('Warning: Mass added')

                Jbar = ((xi(i+1)-xi(i))/2);
                Ueta(i,1) = sum( u.*Gw*Jbar);
                BC = Ueta;
            end
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'one'

        if arg1 || arg2
            Qxi = zeros(1,N);
            for j = 1:N
                Xi  =  ((-1)*arg1+(+1)*arg2)*ones(n,1);
                Eta = (etaAB(j)+etaAB(j+1))/2+(etaAB(j+1)-etaAB(j))/2*etabar;

                [X,Y,dxdxi,dxdeta,dydxi,dydeta] = coordinatemap(Xi,Eta,Domain,DomInfo);

                [qx qy] = exact_solution(X,Y,FunctionType,OutputSelect);

                qxi  = -dxdeta.*qy + dydeta.*qx;

                Jbar = (etaAB(j+1)-etaAB(j))/2;
                Qxi(1,j) = sum( qxi.*Gw*Jbar );
                BC = Qxi;
            end

        elseif arg3 || arg4
            Qeta = zeros(N,1);
            for i = 1:N
                Xi  = (xiLR(i)+xiLR(i+1))/2+(xiLR(i+1)-xiLR(i))/2*xibar;
                Eta =  ((-1)*arg3+(+1)*arg4)*ones(n,1);

                [X,Y,dxdxi,dxdeta,dydxi,dydeta] = coordinatemap(Xi,Eta,Domain,DomInfo);

                [qx qy] = exact_solution(X,Y,FunctionType,OutputSelect);

                qeta =   dxdxi.*qy -  dydxi.*qx;

                Jbar = ((xiLR(i+1)-xiLR(i))/2);
                Qeta(i,1) = sum( qeta.*Gw*Jbar);
                BC = Qeta;
            end
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end