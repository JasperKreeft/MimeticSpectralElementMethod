% function GetLibrary(in)

switch in
    case 'start'
        
        if ispc
            path(path,'G:\MSEM\MSEM_codes\Single_Grid_V5\Library_SingleGrid')
            path(path,'G:\MSEM\MSEM_codes\Single_Grid_V5\HodgeLaplace\Library_HodgeLaplace')
            path(path,'G:\MSEM\MSEM_codes\Single_Grid_V5\HodgeLaplace\Multi_Element\Library_MultiElement')
            path(path,'G:\MSEM\MSEM_codes\Single_Grid_V5\HodgeLaplace\Multi_Element\Darcy\Library_Darcy')
        elseif isunix
            path(path,'/media/My Passport/MSEM/MSEM_codes/Single_Grid_V5/Library_SingleGrid')
            path(path,'/media/My Passport/MSEM/MSEM_codes/Single_Grid_V5/HodgeLaplace/Library_HodgeLaplace')
            path(path,'/media/My Passport/MSEM/MSEM_codes/Single_Grid_V5/HodgeLaplace/Multi_Element/Library_MultiElement')
            path(path,'/media/My Passport/MSEM/MSEM_codes/Single_Grid_V5/HodgeLaplace/Multi_Element/Darcy/Library_Darcy')
        end

    case 'finish'

        if ispc
            rmpath('G:\MSEM\MSEM_codes\Single_Grid_V5\Library_SingleGrid')
            rmpath('G:\MSEM\MSEM_codes\Single_Grid_V5\HodgeLaplace\Library_HodgeLaplace')
            rmpath('G:\MSEM\MSEM_codes\Single_Grid_V5\HodgeLaplace\Multi_Element\Library_MultiElement')
            rmpath('G:\MSEM\MSEM_codes\Single_Grid_V5\HodgeLaplace\Multi_Element\Darcy\Library_Darcy')
        elseif isunix
            rmpath('/media/My Passport/MSEM/MSEM_codes/Single_Grid_V5/Library_SingleGrid')
            rmpath('/media/My Passport/MSEM/MSEM_codes/Single_Grid_V5/HodgeLaplace/Library_HodgeLaplace')
            rmpath('/media/My Passport/MSEM/MSEM_codes/Single_Grid_V5/HodgeLaplace/Multi_Element/Library_MultiElement')
            rmpath('/media/My Passport/MSEM/MSEM_codes/Single_Grid_V5/HodgeLaplace/Multi_Element/Darcy/Library_Darcy')
        end

end