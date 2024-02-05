%CLEARALLBUT  Clear all variables except these.
%   CLEARALLBUT VAR1 removes all variables from the workspace except VAR1.
%
%   CLEARALLBUT VAR1 VAR2 VAR3 ...
%   CLEARALLBUT('VAR1','VAR2','VAR3',...)
%   CLEARALLBUT({'VAR1','VAR2','VAR3',...})
%   do the same for various variable names.
function clearallbut(varargin)
if nargin==0
   error('Expecting an input argument.')
elseif nargin==1
   if ischar(varargin{1})
      in = varargin;
   else
      in = varargin{1};
   end
else
   in = varargin;
end
if ~iscellstr(in)
   error('Input arguments must be strings.')
end
s = evalin('caller','who');
s = setdiff(s,in);
if ~isempty(s)
   %s = sprintf('''%s'',',s{:});
   %s = s(1:end-1);
   %evalin('caller',['clear(',s,')']);
   s = sprintf(' %s',s{:});
   evalin('caller',['clear',s]);
end
