% intobinary    Converts a decimal into a binary vector.
%
% p      decimal representation
% m      length of desired binary vector
%
% v      binary vector
%
% AJT (12/9/18)

function v = intobinary(p,m)

lmax = floor(log2(p));
v = zeros(m,1);
for l = m-1:-1:0
    if (p>=2^l && p<2^(l+1))
        v(m-l,1) = 1;
        p = p - 2^l;
    end
end