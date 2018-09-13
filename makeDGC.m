function K = makeDGC(re,m)

[i j] = upper_indices(m,re);

for p = 1:length(i)
    
    K(:,:,p) = zeros(m,m);
    K(i(p),j(p),p) = 1;
    K(j(p),i(p),p) = 1;
    
end