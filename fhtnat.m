
%------------------------------------------------------
%1D Natural(Hadamard)ordered Fast Hadamard Transform
%------------------------------------------------------
% Author: Gylson Thomas
% e-mail: gylson_thomas@yahoo.com
% Asst. Professor, Electrical and Electronics Engineering Dept.
% MES College of Engineering Kuttippuram,
% Kerala, India, February 2005.
% copyright 2007.

function x=fhtnat(data)
% The function implement the 1D natural(Hadamard)ordered Fast Hadamard Transform,
N = pow2(floor(log2(length(data))));
x = data(1:N);
k1=N; k2=1; k3=N/2;
for i1=1:log2(N)
    L1=1;
    for i2=1:k2
        for i3=1:k3
            i=i3+L1-1; j=i+k3;
            temp1= x(i); temp2 = x(j); 
            x(i) = temp1 + temp2;
            x(j) = temp1 - temp2;
        end
            L1=L1+k1;
    end
        k1 = k1/2;  k2 = k2*2;  k3 = k3/2;
end
x=inv(N)*x; %Delete this line for inverse transform
