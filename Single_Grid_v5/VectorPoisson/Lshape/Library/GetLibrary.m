% function GetLibrary(in)

switch in
    case 'start'
        
        if ispc
            path(path,'C:\Users\Zolder\Documents\MATLAB\MSEM\Single_Grid_Final')
            path(path,'C:\Users\Zolder\Documents\MATLAB\MSEM\Single_Grid_Final\VectorPoisson\Lshape\Library')
        elseif isunix

        end

    case 'finish'

        if ispc
            rmpath('C:\Users\Zolder\Documents\MATLAB\MSEM\Single_Grid_Final')
            rmpath('C:\Users\Zolder\Documents\MATLAB\MSEM\Single_Grid_Final\VectorPoisson\Lshape')
        elseif isunix

        end

end