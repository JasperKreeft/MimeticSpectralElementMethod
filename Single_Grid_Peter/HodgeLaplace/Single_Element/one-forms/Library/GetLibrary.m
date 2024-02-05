% function GetLibrary(in)

switch in
    case 'start'
        
        if ispc
            path(path,'O:\MSEM\MSEM_codes\Single_Grid_Peter\Library_SingleGrid')
            path(path,'O:\MSEM\MSEM_codes\Single_Grid_Peter\HodgeLaplace\Library_HodgeLaplace')
            path(path,'O:\MSEM\MSEM_codes\Single_Grid_Peter\HodgeLaplace\Single_Element\Library_SingleElement')
            path(path,'O:\MSEM\MSEM_codes\Single_Grid_Peter\HodgeLaplace\Single_Element\one-forms\Library')
        elseif isunix
            path(path,'/media/My Passport/MSEM/MSEM_codes/Single_Grid_Peter/Library_SingleGrid')
            path(path,'/media/My Passport/MSEM/MSEM_codes/Single_Grid_Peter/HodgeLaplace/Library_HodgeLaplace')
            path(path,'/media/My Passport/MSEM/MSEM_codes/Single_Grid_Peter/HodgeLaplace/Single_Element/Library_SingleElement')
            path(path,'/media/My Passport/MSEM/MSEM_codes/Single_Grid_Peter/HodgeLaplace/Single_Element/one-forms/Library')
        end

    case 'finish'

        if ispc
            rmpath('O:\MSEM\MSEM_codes\Single_Grid_Peter\Library_SingleGrid')
            rmpath('O:\MSEM\MSEM_codes\Single_Grid_Peter\HodgeLaplace\Library_HodgeLaplace')
            rmpath('O:\MSEM\MSEM_codes\Single_Grid_Peter\HodgeLaplace\Single_Element\Library_SingleElement')
            rmpath('O:\MSEM\MSEM_codes\Single_Grid_Peter\HodgeLaplace\Single_Element\one-forms\Library')
        elseif isunix
            rmpath('/media/My Passport/MSEM/MSEM_codes/Single_Grid_Peter/Library_SingleGrid')
            rmpath('/media/My Passport/MSEM/MSEM_codes/Single_Grid_Peter/HodgeLaplace/Library_HodgeLaplace')
            rmpath('/media/My Passport/MSEM/MSEM_codes/Single_Grid_Peter/HodgeLaplace/Single_Element/Library_SingleElement')
            rmpath('/media/My Passport/MSEM/MSEM_codes/Single_Grid_Peter/HodgeLaplace/Single_Element/one-forms/Library')
        end

end