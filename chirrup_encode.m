%chirrup_encode  Generates K random messages and performs CHIRRUP encoding
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
%
% Y            Y{p} is a 2^m x 2^p matrix of measurements for patch p
% input_bits   B x K matrix of the K B-bit messages
% parity       parity check codes generated in the tree encoding
%
% No. of messages is B = 2^r*[(m-r-p)(m-r-p+3)/2+p-1]-sum(l)  for complex 
%                          B = 2^r*[(m-r-p)(m-r-p+1)/2+p-1)-sum(l) for real
%
% AJT (12/9/18)

function [Y input_bits parity] = chirrup_encode(r,l,re,m,p,K,EbN0)

patches = 2^r;
parity = [];

%calculate length of messages
if (re==0)
    B_patch = m*(m+3)/2 + p - 1;
else
    B_patch = m*(m+1)/2 + p - 1;
end
B = patches*B_patch - sum(l(2:end));

%generate random messages
input_bits = rand(B,K)>0.5;

%tree encoding
if (patches>1)
    [patch_bits parity] = tree_encoder(input_bits,B_patch,patches,l);
    patch_bits = permute(patch_bits,[2 1 3]);
else
    patch_bits = input_bits';
end
    
flag = false;
   
%generate measurements for each patch
for patch = 1:patches

    sigma = sqrt(patches*2^m/(B*EbN0));
    Y{patch} = sim_from_bits(re,m,K,p,sigma,patch_bits(:,:,patch));
    
end