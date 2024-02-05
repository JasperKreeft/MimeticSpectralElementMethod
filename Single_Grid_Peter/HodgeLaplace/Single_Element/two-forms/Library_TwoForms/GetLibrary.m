% function GetLibrary(in)

switch in
    case 'start'
        
        if ispc
            path(path,'O:\msem\MSEM_codes\Single_Grid_Peter\sem')
            path(path,'O:\msem\MSEM_codes\Single_Grid_Peter\Library_SingleGrid')
            path(path,'O:\msem\MSEM_codes\Single_Grid_Peter\HodgeLaplace\Library_HodgeLaplace')
            path(path,'O:\msem\MSEM_codes\Single_Grid_Peter\HodgeLaplace\Single_Element\Library_SingleElement')
            path(path,'O:\msem\MSEM_codes\Single_Grid_Peter\HodgeLaplace\Single_Element\two-forms\Library_TwoForms')
        elseif isunix
            path(path,'/media/My Passport/MSEM/MSEM_codes/Single_Grid_Peter/sem')
            path(path,'/media/My Passport/MSEM/MSEM_codes/Single_Grid_Peter/Library_SingleGrid')
            path(path,'/media/My Passport/MSEM/MSEM_codes/Single_Grid_Peter/HodgeLaplace/Library_HodgeLaplace')
            path(path,'/media/My Passport/MSEM/MSEM_codes/Single_Grid_Peter/HodgeLaplace/Single_Element/Library_SingleElement')
            path(path,'/media/My Passport/MSEM/MSEM_codes/Single_Grid_Peter/HodgeLaplace/Single_Element/two-forms/Library_TwoForms')
        end

    case 'finish'

        if ispc
            rmpath('O:\msem\MSEM_codes\Single_Grid_Peter\sem')
            rmpath('O:\msem\MSEM_codes\Single_Grid_Peter\Library_SingleGrid')
            rmpath('O:\msem\MSEM_codes\Single_Grid_Peter\HodgeLaplace\Library_HodgeLaplace')
            rmpath('O:\msem\MSEM_codes\Single_Grid_Peter\HodgeLaplace\Single_Element\Library_SingleElement')
            rmpath('O:\msem\MSEM_codes\Single_Grid_Peter\HodgeLaplace\Single_Element\two-forms\Library_TwoForms')
        elseif isunix
            rmpath('/media/My Passport/MSEM/MSEM_codes/Single_Grid_Peter/sem')
            rmpath('/media/My Passport/MSEM/MSEM_codes/Single_Grid_Peter/Library_SingleGrid')
            rmpath('/media/My Passport/MSEM/MSEM_codes/Single_Grid_Peter/HodgeLaplace/Library_HodgeLaplace')
            rmpath('/media/My Passport/MSEM/MSEM_codes/Single_Grid_Peter/HodgeLaplace/Single_Element/Library_SingleElement')
            rmpath('/media/My Passport/MSEM/MSEM_codes/Single_Grid_Peter/HodgeLaplace/Single_Element/two-forms/Library_TwoForms')
        end

end