function [P,b] = makePb(re,M,bits)

if (re==0)
    nMuse = M*(M+1)/2;
else
    nMuse = M*(M-1)/2;
end

basis = makeDGC(re,M);

Pbits = bits(1:nMuse);

P = mod( sum(basis(:,:,find(Pbits)),3), 2);

b = bits(nMuse+1:nMuse+M);