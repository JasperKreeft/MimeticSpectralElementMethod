% function GetLibrary(in)

switch in
    case 'start'
        
        if ispc
            disc = pwd;
            path(path,'C:\Users\Jasper\Documents\MATLAB\MSEM_codes\Single_Grid_V5\Library_SingleGrid')
            path(path,'C:\Users\Jasper\Documents\MATLAB\MSEM_codes\Single_Grid_V5\Stokes\Library_Stokes')
            path(path,'C:\Users\Jasper\Documents\MATLAB\MSEM_codes\Single_Grid_V5\Stokes\Multi_Element\Library')
%             disc = pwd;
%             path(path,[disc(1) ':\msem\MSEM_codes\Single_Grid_V5\Library_SingleGrid'])
%             path(path,[disc(1) ':\msem\MSEM_codes\Single_Grid_V5\Stokes\Library_Stokes'])
%             path(path,[disc(1) ':\MSEM\MSEM_codes\Single_Grid_V5\Stokes\Multi_Element\Library'])
        elseif isunix
            path(path,'/media/My Passport/MSEM/MSEM_codes/Single_Grid_V5/Library_SingleGrid')
            path(path,'/media/My Passport/MSEM/MSEM_codes/Single_Grid_V5/Stokes/Library_Stokes')
            path(path,'/media/My Passport/MSEM/MSEM_codes/Single_Grid_V5/Stokes/Multi_Element/Library')
        end

    case 'finish'

        if ispc
            rmpath('C:\Users\Jasper\Documents\MATLAB\MSEM_codes\Single_Grid_V5\Library_SingleGrid')
            rmpath('C:\Users\Jasper\Documents\MATLAB\MSEM_codes\Single_Grid_V5\Stokes\Library_Stokes')
            rmpath('C:\Users\Jasper\Documents\MATLAB\MSEM_codes\Single_Grid_V5\Stokes\Multi_Element\Library')
%             disc = pwd;
%             rmpath([disc(1) ':\msem\MSEM_codes\Single_Grid_V5\Library_SingleGrid'])
%             rmpath([disc(1) ':\msem\MSEM_codes\Single_Grid_V5\Stokes\Library_Stokes'])
%             rmpath([disc(1) ':\MSEM\MSEM_codes\Single_Grid_V5\Stokes\Multi_Element\Library'])
        elseif isunix
            rmpath('/media/My Passport/MSEM/MSEM_codes/Single_Grid_V5/Library_SingleGrid')
            rmpath('/media/My Passport/MSEM/MSEM_codes/Single_Grid_V5/Stokes/Library_Stokes')
            rmpath('/media/My Passport/MSEM/MSEM_codes/Single_Grid_V5/Stokes/Multi_Element/Library')
        end

end