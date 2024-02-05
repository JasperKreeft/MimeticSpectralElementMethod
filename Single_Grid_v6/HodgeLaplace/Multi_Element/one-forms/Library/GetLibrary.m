% function GetLibrary(in)

switch in
    case 'start'

        if ispc
            path(path,'C:\Users\Zolder\Documents\MATLAB\MSEM\Single_Grid_v6\Library_SingleGrid')
            path(path,'C:\Users\Zolder\Documents\MATLAB\MSEM\Single_Grid_v6\HodgeLaplace\Library_HodgeLaplace')
            path(path,'C:\Users\Zolder\Documents\MATLAB\MSEM\Single_Grid_v6\HodgeLaplace\Multi_Element\Library_MultiElement')
            path(path,'C:\Users\Zolder\Documents\MATLAB\MSEM\Single_Grid_v6\HodgeLaplace\Multi_Element\one-forms\Library')
        elseif isunix
            
        elseif ismac
            path(path,'/Users/jasperkreeft/Documents/MATLAB/Codes/MSEM/Single_Grid_v6/Library_SingleGrid')
            path(path,'/Users/jasperkreeft/Documents/MATLAB/Codes/MSEM/Single_Grid_v6/HodgeLaplace/Library_HodgeLaplace')
            path(path,'/Users/jasperkreeft/Documents/MATLAB/Codes/MSEM/Single_Grid_v6/HodgeLaplace/Multi_Element/Library_MultiElement')
            path(path,'/Users/jasperkreeft/Documents/MATLAB/Codes/MSEM/Single_Grid_v6/HodgeLaplace/Multi_Element/one-forms/Library')
        end

    case 'finish'

        if ispc
            rmpath('C:\Users\Zolder\Documents\MATLAB\MSEM\Single_Grid_v6\Library_SingleGrid')
            rmpath('C:\Users\Zolder\Documents\MATLAB\MSEM\Single_Grid_v6\HodgeLaplace\Library_HodgeLaplace')
            rmpath('C:\Users\Zolder\Documents\MATLAB\MSEM\Single_Grid_v6\HodgeLaplace\Multi_Element\Library_MultiElement')
            rmpath('C:\Users\Zolder\Documents\MATLAB\MSEM\Single_Grid_v6\HodgeLaplace\Multi_Element\one-forms\Library')
        elseif isunix
            
        elseif ismac
            rmpath('/Users/jasperkreeft/Documents/MATLAB/Codes/MSEM/Single_Grid_v6/Library_SingleGrid')
            rmpath('/Users/jasperkreeft/Documents/MATLAB/Codes/MSEM/Single_Grid_v6/HodgeLaplace/Library_HodgeLaplace')
            rmpath('/Users/jasperkreeft/Documents/MATLAB/Codes/MSEM/Single_Grid_v6/HodgeLaplace/Multi_Element/Library_MultiElement')
            rmpath('/Users/jasperkreeft/Documents/MATLAB/Codes/MSEM/Single_Grid_v6/HodgeLaplace/Multi_Element/one-forms/Library')
        end

end