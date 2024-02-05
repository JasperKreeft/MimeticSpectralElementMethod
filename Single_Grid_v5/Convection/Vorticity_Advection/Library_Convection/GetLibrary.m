% function GetLibrary(in)

switch in
    case 'start'
        
        if ispc
            path(path,'C:\Users\Jasper\Documents\MATLAB\MSEM_codes\Single_Grid_V5\Library_SingleGrid')
            path(path,'C:\Users\Jasper\Documents\MATLAB\MSEM_codes\Single_Grid_V5\Convection\Vorticity_Advection\Library_Convection')
        elseif isunix
            path(path,'/home/jjkreeft/Convection/Library_SingleGrid')
            path(path,'/home/jjkreeft/Convection/Convection/Vorticity_Advection/Library_Convection')
        end

    case 'finish'

        if ispc
            rmpath('C:\Users\Jasper\Documents\MATLAB\MSEM_codes\Single_Grid_V5\Library_SingleGrid')
            rmpath('C:\Users\Jasper\Documents\MATLAB\MSEM_codes\Single_Grid_V5\Convection\Vorticity_Advection\Library_Convection')
        elseif isunix
            rmpath('/home/jjkreeft/Convection/Library_SingleGrid')
            rmpath('/home/jjkreeft/Convection/Convection/Vorticity_Advection/Library_Convection')
        end

end