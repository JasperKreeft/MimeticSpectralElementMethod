% function GetLibrary(in)

switch in
    case 'start'
        
        if ispc
            path(path,'O:\msem\MSEM_codes\Single_Grid\Library_SingleGrid')
            path(path,'O:\MSEM\MSEM_codes\Single_Grid\VectorPoisson\Lshape\Library')
        elseif isunix
            path(path,'/media/FREECOM HDD/msem/MSEM_codes/Single_Grid/Library_SingleGrid')
            path(path,'/media/FREECOM HDD/msem/MSEM_codes/Single_Grid/VectorPoisson/Lshape/Library')
        end

    case 'finish'

        if ispc
            rmpath('O:\msem\MSEM_codes\Single_Grid\Library_SingleGrid')
            rmpath('O:\msem\MSEM_codes\Single_Grid\VectorPoisson\Lshape\Library')
        elseif isunix
            rmpath('/media/FREECOM HDD/msem/MSEM_codes/Single_Grid/Library_SingleGrid')
            rmpath('/media/FREECOM HDD/msem/MSEM_codes/Single_Grid/VectorPoisson/Lshape/Library')
        end

end