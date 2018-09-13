% chirrup_test  Runs multiple trials of CHIRRUP for given parameters.
%
% r          0, 1 or 2; 2^r patches
% l          length-2^r vector of parity check digits
%            recommend: r=0: l=0, r=1: l=[0 15], r=2: l = [0 10 10 15]
% re         logical: false=complex chirps, true=real chirps
% m          size of chirp code (2^m is length of codewords in each slot)
%            recommend m = 6, 7 or 8
% p          2^p slots
%            require p<=(m-r)(m-r+3)/2-1 for complex
%                    p<=(m-r)(m-r+1)/2-1 for real         
% K          number of messages
% EbN0       energy-per-bit (Eb/N0)
% trials     number of trials
%
% propfound    proportion of messages recovered
% ave_time     average running time per reconstruction
%
% No. of messages is B = 2^r*[(m-r-p)(m-r-p+3)/2+p-1]-sum(l)  for complex 
%                          B = 2^r*[(m-r-p)(m-r-p+1)/2+p-1)-sum(l) for real
%
% AJT (12/9/18)

function [propfound ave_time] = chirrup_test(r,l,re,m,p,K,EbN0,trials)

sumpropfound = 0;
sumtiming = 0;

for t = 1:trials
    
    [Y input_bits parity] = chirrup_encode(r,l,re,m,p,K,EbN0);
    
    [output_bits timing_trial] = chirrup_decode(Y,r,l,parity,re,m,p,K);
    
    propfound_trial = compare_bits(input_bits,output_bits);
    
    sumpropfound = sumpropfound + propfound_trial;
    sumtiming = sumtiming + timing_trial;
    
end

propfound = sumpropfound/trials;
ave_time = sumtiming/trials;