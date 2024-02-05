% function GetLibrary(in)

switch in
    case 'start'
        
        if ispc
            path(path,'C:\Users\Zolder\Documents\MATLAB\MSEM\Single_Grid_v6\Library_SingleGrid')
            path(path,'C:\Users\Zolder\Documents\MATLAB\MSEM\Single_Grid_v6\Stokes\Library_Stokes')
            path(path,'C:\Users\Zolder\Documents\MATLAB\MSEM\Single_Grid_v6\Stokes\Multi_Element\Library')
        elseif isunix

        end

    case 'finish'

        if ispc
            rmpath('C:\Users\Zolder\Documents\MATLAB\MSEM\Single_Grid_v6\Library_SingleGrid')
            rmpath('C:\Users\Zolder\Documents\MATLAB\MSEM\Single_Grid_v6\Stokes\Library_Stokes')
            rmpath('C:\Users\Zolder\Documents\MATLAB\MSEM\Single_Grid_v6\Stokes\Multi_Element\Library')
        elseif isunix

        end

end