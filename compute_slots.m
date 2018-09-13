function comps = compute_slots(bits,re,M,p)
if (re==0)
    nMuse = M*(M+1)/2;
else
    nMuse = M*(M-1)/2;
end
comps(1) = outofbinary(bits(end-p+1:end))+1;

trans = bits(nMuse+M-1:-1:nMuse+M-p);
if outofbinary(trans)==0
    trans(1) = 1;
end
comps(2) = outofbinary(mod(bits(end-p+1:end)+trans,2))+1;