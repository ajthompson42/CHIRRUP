% outofbinary    Converts a binary vector into decimal.
%
% v      binary vector
%
% m      decimal representation
%
% AJT (12/9/18)

function m = outofbinary(v)

v = flipud(v);
m = 0;
mult = 1;
for l = length(v):-1:1
    m = m + v(l)*mult;
    mult = mult*2;
end