function [ out1 out2 ] = exact_solution(X,Y,FunctionType,OutputSelect)

global Re

%% Exact Solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch FunctionType
    
    case 'reguralizedLDC'
        
        [psi,u,v,w,p,fy]=ReguralizedLDC(Re,X,Y);
        
        switch OutputSelect
            case 'streamfunction'
                out1 = psi;
            case 'vorticity'
                out1 = w;
            case 'velocity'
                out1 = u;
                out2 = v;
            case 'pressure'
                out1 = p;
            case 'force'
                fx = 0*X;
                out1 = fx;
                out2 = fy;
            case 'mass'
                out1 = 0;
        end

    case 'reguralizedLDC_p'
        [Fx,Fy]=ReguralizedLDC_p(Re,X,Y);
                out1 = Fx;
                out2 = Fy;
end