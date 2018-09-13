# CHIRRUP
CHIRRUP algorithm for unsourced multiple access

This toolbox performs encoding and decoding using CHIRRUP for unsourced 
multiple access, as described in the paper 'CHIRRUP: a practical algorithm 
for unsourced multiple access' (Calderbank/Thompson 2018). There is also 
code for performing repeated tests to obtain average success proportions 
and timings.

The main programs are:
chirrup_encode: performs CHIRRUP encoding of random messages
chirrup_decode: performs CHIRRUP decoding
chirrup_test: performs repeated tests

To perform CHIRRUP encoding,
[Y,input_bits] = chirrup_encode(r,l,re,m,p,K,EbN0)

For example,
[Y input_bits] = chirrup_encode(0,0,1,8,7,200,50)
performs CHIRRUP encoding of 200 messages as 1=2^0 patches (with no parity
bits trivially), using real binary chirps of size 2^8 in 2^7 slots, with
energy-per-bit in Eb/N0=50.

[Y input_bits parity] = chirrup_encode(1,[0 15],0,7,7,100,20)
performs CHIRRUP encoding of 100 messages as 2=2^1 patches (with 15 parity
bits), using complex binary chirps of size 2^7 in 2^7 slots, with
energy-per-bit Eb/N0=20.

To perform CHIRRUP decoding,
[output_bits timing] = chirrup_decode(Y,r,l,parity,re,m,p,K,params)

For example,
[output_bits timing] = chirrup_decode(Y,0,0,[],1,8,7,200)
performs CHIRRUP decoding of 200 messages as 1=2^0 patches (with no parity
bits trivially), using real binary chirps of size 2^8 in 2^7 slots.

[output_bits timing] = chirrup_decode(Y,1,[0 15],parity,0,7,7,100)
performs CHIRRUP decoding of 100 messages as 2=2^1 patches (with 15 parity
bits), using complex binary chirps of size 2^7 in 2^7 slots.

The proportion of correctly identified messages can be found using
propfound = compare_bits(input_bits,output_bits).

A repeated test over a number of trials can be run using
[propfound ave_time] = chirrup_test(r,l,re,m,p,K,EbN0,trials).

The code is largely authored by Andrew Thompson (National Physical 
Laboratory, UK), but some of the original chirp reconstruction code is also 
written by Sina Jafarpour (Facebook, USA). The fast Hadamard routine is by
Gylson Thomas (Jyothi Engineering College, India).

Please refer to the code for more detailed explanation of input and output 
parameters for each program.
