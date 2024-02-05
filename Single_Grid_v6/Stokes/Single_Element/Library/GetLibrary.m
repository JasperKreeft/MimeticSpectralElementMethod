% function GetLibrary(in)

switch in
    case 'start'
        
        if ispc
            path(path,'C:\Users\jaspe\Documents\MATLAB\Codes\MSEM\Single_Grid_v6\Library_SingleGrid')
            path(path,'C:\Users\jaspe\Documents\MATLAB\Codes\MSEM\Single_Grid_v6\Stokes\Library_Stokes')
            path(path,'C:\Users\jaspe\Documents\MATLAB\Codes\MSEM\Single_Grid_v6\Stokes\Single_Element\Library')
        elseif isunix

        end

    case 'finish'

        if ispc
            rmpath('C:\Users\jaspe\Documents\MATLAB\Codes\MSEM\Single_Grid_v6\Library_SingleGrid')
            rmpath('C:\Users\jaspe\Documents\MATLAB\Codes\MSEM\Single_Grid_v6\Stokes\Library_Stokes')
            rmpath('C:\Users\jaspe\Documents\MATLAB\Codes\MSEM\Single_Grid_v6\Stokes\Single_Element\Library')
        elseif isunix

        end

end