%{
DESCRIPTION:
A convenience alias for tensor

SYNTAX:
Same as tensor
 
AUTHOR: 
Mario Merino <mario.merino@uc3m.es>
%}
function T = scalar(varargin)
    T = anakin.tensor(varargin{:});
end