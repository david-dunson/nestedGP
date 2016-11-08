function model = mysetparam(model, psi)

%@SSMODEL/SETPARAM Set model parameter values.
%   model = SETPARAM(model, psi[, transformed])

% (c) 2006-2007 Jyh-Ying Peng ´^´¼·­
% $Revision 1.0.0 $ $Generated: 2007/09/04 $


    model.psi       = set(model.psi, psi);
    psi             = model.psi.value;



