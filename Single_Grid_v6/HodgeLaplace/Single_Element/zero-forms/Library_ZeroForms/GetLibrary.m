% function GetLibrary(in)

switch in
    case 'start'

        if ispc
            path(path,'..\..\..\Library_SingleGrid')
            path(path,'..\..\Library_HodgeLaplace')
            path(path,'..\Library_SingleElement')
            path(path,'Library_ZeroForms')
        elseif isunix

        end

    case 'finish'

        if ispc
            rmpath('..\..\..\Library_SingleGrid')
            rmpath('..\..\Library_HodgeLaplace')
            rmpath('..\Library_SingleElement')
            rmpath('Library_ZeroForms')
        elseif isunix

        end

end