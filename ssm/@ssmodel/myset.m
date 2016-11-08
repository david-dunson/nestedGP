function model = myset(model, A, M)

%@SSMODEL/SET Set model elements.
%   model = SET(model0, A, M[, func, grad, psi, pmask])

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

if ischar(A), A = {A}; elseif ~iscell(A), error('ssm:ssmodel:set:InputError', 'A must be a string or a cell array of strings.'); end
%if ~isa(M, 'ssmat'), if isnumeric(M), M = ssmat(M); else error('ssm:ssmodel:set:InputError', 'M must be a SSMAT or a matrix.'); end, end


switch A{1}(1)
    case 'H'
        model.H = M;
    case 'Q'
        model.Q = M;
end
       



