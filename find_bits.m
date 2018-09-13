% check_result      compare result of findPbc with simulation
%function OK = check_result(res,sim)

% SJS 14/11/07

function bits = find_bits(res,re,M,p)

bits = [];
OK = [];
nOK = 0;

K = length(res);
s = 1;

for k = 1:K
    
    Pk = res(k).P;
    bk = res(k).b;
    ck = res(k).c;
    compsk = res(k).comps;
    [i j] = upper_indices(M,re);
    for l = 1:length(i)
        Pvec(l) = Pk(i(l),j(l));
    end
    bits2 = [Pvec bk];
    bits2 = bits2(1,2:end);
    slot_bit = (intobinary(compsk(1)-1,p))';
    bits(k,:) = [bits2 slot_bit];
end