% genK      generate a binary chirp from P & b
%           
% Sina Jafarpour (28/3/08) and AJT (13/9/18)

function rm = gen_chirp(P,b)

M = length(b);

% construct Reed-Muller code from P and b
rm = zeros(2^M,1);
a = zeros(M,1);
for q = 1:2^M
    sum1 = a'*P*a;  
    sum2 = b*a;
    rm(q) = i^sum1 * (-1)^sum2;
    % next a
    for ix = M:-1:1
        if a(ix)==1
            a(ix)=0;
        else
            a(ix)=1;
            break;
        end
    end
end