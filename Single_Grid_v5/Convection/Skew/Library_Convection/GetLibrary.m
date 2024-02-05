% function GetLibrary(in)

switch in
    case 'start'
        
        if ispc
            path(path,'C:\Users\Jasper\Documents\MATLAB\MSEM_codes\Single_Grid_V5\Library_SingleGrid')
            path(path,'C:\Users\Jasper\Documents\MATLAB\MSEM_codes\Single_Grid_V5\Convection\Skew\Library_Convection')
        elseif isunix
            path(path,'/media/My Passport/MSEM/MSEM_codes/Single_Grid_V5/Library_SingleGrid')
            path(path,'/media/My Passport/MSEM/MSEM_codes/Single_Grid_V5/Convection/Skew/Library_Convection')
        end

    case 'finish'

        if ispc
            rmpath('C:\Users\Jasper\Documents\MATLAB\MSEM_codes\Single_Grid_V5\Library_SingleGrid')
            rmpath('C:\Users\Jasper\Documents\MATLAB\MSEM_codes\Single_Grid_V5\Convection\Skew\Library_Convection')
        elseif isunix
            rmpath('/media/My Passport/MSEM/MSEM_codes/Single_Grid_V5/Library_SingleGrid')
            rmpath('/media/My Passport/MSEM/MSEM_codes/Single_Grid_V5/Convection/Skew/Library_Convection')
        end

end