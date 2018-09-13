%chirrup_decode  Performs CHIRRUP decoding
%
% Y          Y{p} is a 2^m x 2^p matrix of measurements for patch p
% r          0, 1 or 2; 2^r patches
% l          length-2^r vector of parity check digits
%            recommend: r=0: l=0, r=1: l=[0 15], r=2: l = [0 10 10 15]
% parity     parity check codes generated in the tree encoding
%            leave empty if only 1 patch (r=0) is used
% re         logical: false=complex chirps, true=real chirps
% m          size of chirp code (2^m is length of codewords in each slot)
%            recommend m = 6, 7 or 8
% p          2^p slots
%            require p<=(m-r)(m-r+3)/2-1 for complex
%                    p<=(m-r)(m-r+1)/2-1 for real         
% K          number of messages
% params     various; see below
%
% output_bits  B x K matrix of the K B-bit messages
% timing       running time in seconds
%
% AJT (12/9/18)

function [output_bits timing] = chirrup_decode(Y,r,l,parity,re,m,p,K,params_in)

%default parameter values
params.alpha = 0.1; %accept components if coefficient is within alpha of 1
params.circuits = 5; %number of circuits in peeling decoder
params.sparsity_factor = 3; %factor by which number of iterations exceeds 
%expected sparsity for a given slot
params.tree_order = 3; %number of candidates per row in findPb tree search

%read in parameter inputs
if (nargin==9)
    if isfield(params_in,'alpha')
        params.alpha = params_in.alpha;
    end
    if isfield(params_in,'circuits')
        params.circuits = params_in.circuits;
    end
    if isfield(params_in,'sparsity_ratio')
        params.sparsity_factor = params_in.sparsity_factor;
    end
    if isfield(params_in,'tree_order')
        params.tree_order = params_in.tree_order;
    end 
end

global outer_recov

warning('off','MATLAB:rankDeficientMatrix');

patches = 2^r;
flag = false;

for patch = 1:patches
        
        tic
        outer_recov = [];
        count = 1; 
        Yp = Y{patch};
        found = [];

        %cycle through the slots repeatedly
        for c = 1:params.circuits

            for slot = 1:2^p
                
                %run chirp reconstruction on given slot
                sparsity_ratio = 3;
                recov = chirp_rec(Yp(:,slot),re,ceil(params.sparsity_factor*K/2^(p-1)),slot,params);
                
                [i j] = upper_indices(m,re);

                if (~isempty(recov(1).P))
                    
                    %find slot pairs for each recovered component
                    for r = 1:size(recov,2)
                        Q=recov(r).P;
                        for k = 1:length(i)
                            Pvec(k) = Q(i(k),j(k));
                        end
                        Pbvec = [Pvec recov(r).b];
                        already = 0;
                        if (re==0)
                            trans = Pbvec(m*(m+1)/2+m:-1:m*(m+1)/2+m-p+1);
                        else
                            trans = Pbvec(m*(m-1)/2+m:-1:m*(m-1)/2+m-p+1);
                        end
                        if outofbinary(trans)==0
                            trans(1) = 1;
                        end
                        
                        if (re==0)
                            if recov(r).P(1,1)==1
                                recov(r).P(1,1) = 0;
                                recov(r).comps(2) = slot;
                                recov(r).comps(1) = outofbinary(mod(trans+(intobinary(slot-1,p))',2))+1;
                            else
                                recov(r).comps(1) = slot;
                                recov(r).comps(2) = outofbinary(mod(trans+(intobinary(slot-1,p))',2))+1;
                            end
                        else
                            if recov(r).P(1,2)==1
                                recov(r).P(1,2) = 0;
                                recov(r).P(2,1) = 0;
                                recov(r).comps(2) = slot;
                                recov(r).comps(1) = outofbinary(mod(trans+(intobinary(slot-1,p))',2))+1;
                            else
                                recov(r).comps(1) = slot;
                                recov(r).comps(2) = outofbinary(mod(trans+(intobinary(slot-1,p))',2))+1;
                            end
                        end
                        
                        %accept component if its coefficient is close to 1
                        %if (abs(recov(r).c - 1)<params.alpha) %alternative condition
                        if (real(recov(r).c)>1-params.alpha && real(recov(r).c)<1+params.alpha && abs(imag(recov(r).c))<params.alpha)
                            for r2 = 1:count-1
                                if (min(min(recov(r).P==outer_recov(r2).P)) && min(recov(r).b==outer_recov(r2).b))
                                    already = 1;
                                end
                            end
                            if already==0
                                outer_recov(count).P = recov(r).P;
                                outer_recov(count).b = recov(r).b;
                                outer_recov(count).c = recov(r).c;
                                outer_recov(count).comps = recov(r).comps;
                                count = count + 1;
                            end
                        end
                    end
                end
            end
        end
        
        %convert components back into bits
        bit_matrix = find_bits(outer_recov,re,m,p);
        if (size(bit_matrix,1)==0) 
            flag = 1;
        end
        processed_bits{patch} = bit_matrix';
        timings(patch) = toc;

end
    
    %patch together messages if necessary
    tic
    if (flag==1)
        propfound_trial = 0;
        sumpropfound = sumpropfound + propfound_trial;
        timings(patches+1) = toc;
    else
        if (patches>1)
            output_bits = tree_decoder(processed_bits,l,parity);
        else
            output_bits = bit_matrix';
        end    
        kmin = min(K,size(output_bits,2));
        output_bits = output_bits(:,1:kmin);
        timings(patches+1) = toc;
    end
    
    timing = sum(timings);