% function GetLibrary(in)

switch in
    case 'start'
        
        if ispc
            path(path,'O:\msem\MSEM_codes\Single_Grid_V5\Library_SingleGrid')
            path(path,'O:\msem\MSEM_codes\Single_Grid_V5\Stokes\Library_Stokes')
            path(path,'O:\MSEM\MSEM_codes\Single_Grid_V5\Stokes\Single_Element\Library')
        elseif isunix
            path(path,'/media/My Passport/MSEM/MSEM_codes/Single_Grid_V5/Library_SingleGrid')
            path(path,'/media/My Passport/MSEM/MSEM_codes/Single_Grid_V5/Stokes/Library_Stokes')
            path(path,'/media/My Passport/MSEM/MSEM_codes/Single_Grid_V5/Stokes/Single_Element/Library')
        end

    case 'finish'

        if ispc
            rmpath('O:\msem\MSEM_codes\Single_Grid_V5\Library_SingleGrid')
            rmpath('O:\msem\MSEM_codes\Single_Grid_V5\Stokes\Library_Stokes')
            rmpath('O:\msem\MSEM_codes\Single_Grid_V5\Stokes\Single_Element\Library')
        elseif isunix
            rmpath('/media/My Passport/MSEM/MSEM_codes/Single_Grid_V5/Library_SingleGrid')
            rmpath('/media/My Passport/MSEM/MSEM_codes/Single_Grid_V5/Stokes/Library_Stokes')
            rmpath('/media/My Passport/MSEM/MSEM_codes/Single_Grid_V5/Stokes/Single_Element/Library')
        end

end