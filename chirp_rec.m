% chirp_rec     Runs the chirp reconstruction algorithm, incorporating 
%               components found in other slots.
%
% y        measurement vector for given slot
% re       logical: false=complex chirps, true=real chirps
% nitlim   iteration limit
% slot     slot number
%
% recov    multidimensional struct storing (P,b,c) for each component
%          e.g. for component i: recov(1).P, recov(1).b, recov(1).c
%
% Sina Jafarpour (14/4/08) and AJT (12/9/18)

function recov = chirp_rec(y,re,nitlim,slot,params)

global outer_recov

foo.P = [];
foo.b = [];
foo.c = [];
foo.comps = [];
recov(1) = foo;
allRM = [];
ncomp = 0;

%include components found in other slots
count = 1;
for r = 1:size(outer_recov,2)
    if outer_recov(r).comps(1)==slot
        recov(count).P = outer_recov(r).P;
        recov(count).b = outer_recov(r).b;
        recov(count).c = outer_recov(r).c;
        count = count + 1;
    elseif outer_recov(r).comps(2)==slot
        recov(count).P = outer_recov(r).P;
        if (re==0)
            recov(count).P(1,1) = 1;
        else
            recov(count).P(1,2) = 1;
            recov(count).P(2,1) = 1;
        end
        recov(count).b = outer_recov(r).b;
        recov(count).c = outer_recov(r).c;
        count = count + 1;
    end
end

%peel off components and recalculate residual
if (count>1)
for r = 1:size(recov,2)
    RM = gen_chirp(recov(r).P,recov(r).b); 
    allRM = [allRM RM];
end
cr = allRM\y;
y = y - allRM*cr;
ncomp = size(recov,2);
end

M = log2(length(y));  

while (ncomp < nitlim && norm(y)>1e-3)
    
    [Phat bhat] = findPb(y,re,params);
    
    %determine component from P,b
    RM = gen_chirp(Phat,bhat);  
    allRM = [allRM RM];
    
    %use all components & refine previous ests.
    cr = allRM\y;
    c = cr(end);
    allc = [recov(1:ncomp).c].';
    newc = [allc; 0]+cr;
    for q = 1:ncomp
        recov(q).c = newc(q);
    end
    y = y - allRM*cr;

    %a component has been detected and handled.
    %record component parameters
    ncomp = ncomp+1;
    recov(ncomp).P = Phat;
    recov(ncomp).b = bhat;
    recov(ncomp).c = c;
    recov(ncomp).comp = RM;
    
end