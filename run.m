
function [propfound, ave_time] = run(r, l, re, m, p, EbN0, K, trials)
    addpath('utils');
%master testing function, runs multiple trials, for definitions of
%parameters see the Encoder and Decoder classes

    params_in=[];
    sumpropfound=0;
    sumtiming=0;

    
    for trial = 1:trials
        [propfound_trial, timing_trial] = chirrup_test(r,l,re,m,p,K,EbN0,[],params_in);
        sumpropfound=sumpropfound+propfound_trial;
        sumtiming=sumtiming+timing_trial;
    end

    propfound = sumpropfound/trials;
    ave_time = sumtiming/trials;

end



function [propfound_trial, timing_trial] = chirrup_test(r,l,re,m,p,K,EbN0,input_bits,params_in)
    %runs a single test
    encoder=Encoder(r,l,re,m,p,K,EbN0,input_bits);
    [encoder, input_bits] = encoder.generate_random_bits(); %we update the value of input_bits to save confusion, despite already being encoder.input_bits
    [Y, parity] = encoder.chirrup_encode;
    decoder=Decoder(Y,r,l,parity,re,m,p,K,params_in);
    [output_bits, timing_trial] = decoder.chirrup_decode();
    propfound_trial = compare_bits(input_bits,output_bits);
end
