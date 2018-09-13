% probe_ye  probe data vector with error vector
%
% [prb,yprod] = probe_ye(y,e)
%
% given data vector y and "error value" e,
% compute Hadamard transform prb  of yprod = conj(y(a)) * y(a+e)

% Sina Jafarpour 28/3/08

function [prb,yprod] = probe_ye(y,ein)

Mpow = length(y);
M = log2(Mpow);  % this should be an integer
if isscalar(ein)
    e = ein;
else
    e = (2.^(M-1:-1:0))*ein(:);
end

% index vector (starting at 0)
a = (0:Mpow-1)';

% a plus e values
ape = bitxor(a,e);

% y(a+e)
yape = y(ape+1);

% product
yprod = conj(y).*yape; 

% use fast Walsh-Hadamard transform routine
prb = fhtnat(yprod);