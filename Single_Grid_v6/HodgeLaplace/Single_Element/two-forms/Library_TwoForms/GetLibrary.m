% function GetLibrary(in)

switch in
    case 'start'
        
        if ispc
            path(path,'C:\Users\Zolder\Documents\MATLAB\MSEM\Single_Grid_v6\Library_SingleGrid')
            path(path,'C:\Users\Zolder\Documents\MATLAB\MSEM\Single_Grid_v6\HodgeLaplace\Library_HodgeLaplace')
            path(path,'C:\Users\Zolder\Documents\MATLAB\MSEM\Single_Grid_v6\HodgeLaplace\Single_Element\Library_SingleElement')
            path(path,'C:\Users\Zolder\Documents\MATLAB\MSEM\Single_Grid_v6\HodgeLaplace\Single_Element\two-forms\Library_TwoForms')
        elseif isunix

        end

    case 'finish'

        if ispc
            rmpath('C:\Users\Zolder\Documents\MATLAB\MSEM\Single_Grid_v6\Library_SingleGrid')
            rmpath('C:\Users\Zolder\Documents\MATLAB\MSEM\Single_Grid_v6\HodgeLaplace\Library_HodgeLaplace')
            rmpath('C:\Users\Zolder\Documents\MATLAB\MSEM\Single_Grid_v6\HodgeLaplace\Single_Element\Library_SingleElement')
            rmpath('C:\Users\Zolder\Documents\MATLAB\MSEM\Single_Grid_v6\HodgeLaplace\Single_Element\two-forms\Library_TwoForms')
        elseif isunix

        end

end