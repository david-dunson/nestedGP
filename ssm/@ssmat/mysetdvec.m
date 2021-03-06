function m = mysetdvec(m, dsubvec)

%@SSMAT/SETDVEC Update state space matrix dynamic part.
%   m = SETDVEC(m, dsubvec)

% (c) 2006-2007 Jyh-Ying Peng �^����
% $Revision 1.0.0 $ $Generated: 2007/09/04 $

%%%%%%%% INFO: Update to dvec, if it's previously shorter than n it will
%%%%%%%% expand, with zeros at places not updated at new time points,
%%%%%%%% unchanged elsewhere. If it's previously longer than n, dvec beyond
%%%%%%%% time point n will be unchanged, but also unused since algorithms
%%%%%%%% never reference beyond n.


 m.dvec = dsubvec;
