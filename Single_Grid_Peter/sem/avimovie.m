function []=avimovie(filename,figurename,first_condition,last_condition)
% avimovie is a script to create a avi-file by taking snap-shots of you
% figure. Because avi-files created in Matlab can become really big, the
% mencoder is used, which reduces the file size considerably.
%
% first_condition is 1 when it is the first time you enter this function
% and 0 elsewhere.
% last_condition is 1 when it is the last time you enter this function
% and 0 elsewhere
%
% Written by Jasper Kreeft 2011
%

global aviobj

filename = strcat(filename,'.avi');

if first_condition
    aviobj = avifile(filename,'fps',8,'quality',100,'compression','None');
end

if size(get(gcf,'visible'),2)==2
    frame = getframe(figurename);
    aviobj = addframe(aviobj,frame);
else    
    aviobj = addframe(aviobj,figurename);
end



if last_condition

    close(figurename);
    aviobj = close(aviobj);

    if isunix
    avifilename      = filename;
    avifilename_comp = strcat(filename(1:end-4),'_comp','.avi');

    system(['mencoder ' avifilename ' -o ' avifilename_comp ' -ovc lavc']);
    system(['mv ' avifilename_comp ' ' avifilename]);
    system(['mplayer -fs -fps 8 ' avifilename]);
    end

end