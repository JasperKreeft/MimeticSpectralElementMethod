% function GetLibrary(in)

switch in
    case 'start'
        
        if ispc
%             disk = pwd; disk = disk(1);
%             path(path,[disk ':\msem\MSEM_codes\Single_Grid_V5\Library_SingleGrid'])
%             path(path,[disk ':\msem\MSEM_codes\Single_Grid_V5\Stokes\Library_Stokes'])
%             path(path,[disk
%             ':\MSEM\MSEM_codes\Single_Grid_V5\Stokes\CoVectorValued_multi\Library'])
            path(path,'C:\Users\Jasper\Documents\MATLAB\MSEM_codes\Single_Grid_V5\Library_SingleGrid')
            path(path,'C:\Users\Jasper\Documents\MATLAB\MSEM_codes\Single_Grid_V5\Stokes\Library_Stokes')
            path(path,'C:\Users\Jasper\Documents\MATLAB\MSEM_codes\Single_Grid_V5\Stokes\CoVectorValued_multi\Library')

        elseif isunix
            path(path,'/media/My Passport/MSEM/MSEM_codes/Single_Grid_V5/Library_SingleGrid')
            path(path,'/media/My Passport/MSEM/MSEM_codes/Single_Grid_V5/Stokes/Library_Stokes')
            path(path,'/media/My Passport/MSEM/MSEM_codes/Single_Grid_V5/Stokes/CoVectorValued_multi/Library')
        end

    case 'finish'

        if ispc
%             disk = pwd; disk = disk(1);
%             rmpath([disk ':\msem\MSEM_codes\Single_Grid_V5\Library_SingleGrid'])
%             rmpath([disk ':\msem\MSEM_codes\Single_Grid_V5\Stokes\Library_Stokes'])
%             rmpath([disk ':\msem\MSEM_codes\Single_Grid_V5\Stokes\CoVectorValued_multi\Library'])
            rmpath('C:\Users\Jasper\Documents\MATLAB\MSEM_codes\Single_Grid_V5\Library_SingleGrid')
            rmpath('C:\Users\Jasper\Documents\MATLAB\MSEM_codes\Single_Grid_V5\Stokes\Library_Stokes')
            rmpath('C:\Users\Jasper\Documents\MATLAB\MSEM_codes\Single_Grid_V5\Stokes\CoVectorValued_multi\Library')

        elseif isunix
            rmpath('/media/My Passport/MSEM/MSEM_codes/Single_Grid_V5/Library_SingleGrid')
            rmpath('/media/My Passport/MSEM/MSEM_codes/Single_Grid_V5/Stokes/Library_Stokes')
            rmpath('/media/My Passport/MSEM/MSEM_codes/Single_Grid_V5/Stokes/CoVectorValued_multi/Library')
        end

end