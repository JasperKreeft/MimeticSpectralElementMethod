function f = forcefunction_assembly(form,Mesh,FunctionType)

global numElements
global nr_0 globalnr_0

disp('assembly force function')

switch form
    
    case 0
        
        f  = zeros(nr_0,1);

        for i=1:numElements

        ind0 = globalnr_0(:,i);

        F = forcefunction(0,Mesh.X(:,i),Mesh.Y(:,i),FunctionType);
        f(ind0) = F;

        end
        
        
    case 1
        
        
        
    case 2
        
        
end