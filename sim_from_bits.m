% sim_from_bits  Generates binary chirp measurements from bits
%
% re         logical: false=complex chirps, true=real chirps
% m          size of chirp code (2^m is length of codewords in each slot)
%            recommend m = 6, 7 or 8
% K          number of messages
% p          2^p slots
%            require p<=(m-r)(m-r+3)/2-1 for complex
%                    p<=(m-r)(m-r+1)/2-1 for real  
% sigma      SD of noise: sigma = sqrt(patches*2^m/(B*EbN0))
% bits       k x 2^m matrix of bits to encode
%
% Y          length-2^m vector of measurements
%
% AJT (12/9/18)

function Y = sim_from_bits(re,m,K,p,sigma,bits)

Y = zeros(2^m,2^p);

for k = 1:K
    
    %compute slots
    comps = compute_slots(bits(k,:),re,m,p);
   
    %make (P,b) for each slot
    bits1 = [0 bits(k,:)];
    bits2 = [1 bits(k,:)];
    [Pee1,bee1] = makePb(re,m,bits1);
    [Pee2,bee2] = makePb(re,m,bits2);
    
    %generate binary chirp vector for each slot
    rm1 = gen_chirp(Pee1,bee1);
    rm2 = gen_chirp(Pee2,bee2);
    
    %add onto measurement
    Y(:,comps(1)) = Y(:,comps(1))+rm1;
    Y(:,comps(2)) = Y(:,comps(2))+rm2;
end

%add noise (Gaussian for real, Complex Gaussian for complex)
if (re==0)
    Y = Y + repmat(sigma*(randn(2^m,1)+1i*randn(2^m,1)),[1 2^p]);
else
    Y = Y + repmat(sigma*randn(2^m,1),[1 2^p]);
end