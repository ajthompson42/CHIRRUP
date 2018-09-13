% findPb  Find P matrix & b vector recursively.
%
% y         signal, length 2^M
% re        logical: false=complex chirps, true=real chirps
%
% Pout      P matrix
% bout      b vector
%     
% Sina Jafarpour (29/4/08) and AJT (12/9/18)

function [Pout bout] = findPb(y,re,params)

global perm prout DGbasis;  
% so these are visible from within recursive function
twoM = length(y);
M = log2(twoM);

DGbasis = makeDGC(re,M);

if (re==0)
    DGblen = M*(M+1)/2;
else
    DGblen = M*(M-1)/2;
end

%determine probe response to each basis vector
%probe = Hadamard xform of y^*(a)y(a+e), e = basis vector
evp = eye(M);
prout = zeros(M,twoM);
for m = 1:M
    prout(m,:) = probe_ye(y,evp(m,:));
end

perm = (1:M)';
    
%manipulate permutation so only rows encompassing the first
%DGblen elements will be used.
if (re==0)
    csum = cumsum(M:-1:1);
else
    csum = cumsum(M-1:-1:1);
end
nrow = find(csum>=DGblen);
nrow = nrow(1);
p1 = perm(perm<=nrow);
p2 = perm(perm>nrow);
perm = [p1; p2];   %could have just p1, since we'll never loop over rest

%do the search
[Pout,bout,rec] = fPN_recursive(re,1,y,zeros(M),[],params);

%what to do if nothing found?
if isempty(Pout)
    % Oh heck.  Didn't find anything.
    % Make one up at random...  hopefully this will kick the signal enough
    % to make something stand out next time, and this bogus component will
    % be subsequently removed.
    [Pout,bout] = bstr2Pb(M,[]);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Pout,bout,crec] = fPN_recursive(re,m,y,P,rec,params)

global perm prout DGbasis;  
%defined prior to call of this recursive function

M = length(P);

%which basis matrices go with current row
DGblen = size(DGbasis,3);
if (re==0)
    csM1 = cumsum(M:-1:1);
else
    csM1 = cumsum(M-1:-1:1);
end
DGto = csM1(m);
DGto = min(DGto,DGblen);

pm = perm(m);
row = P(pm,:);
if (re==0)
    fixed = perm(1:m-1);  
    unfixed = perm(m:end);
else
    fixed = perm(1:m);
    unfixed = perm(m+1:end);
end

%generate all candidates for this row,
%given constraints by already-fixed rows
%(= fixed element values in this row)
nfree = 2^length(unfixed); 
                          
nums = 0:nfree-1;
cands = row(ones(nfree,1),:);
cands(:,unfixed) = (dec2bin(nums)=='1');
candix = cands * (2.^(M-1:-1:0)')  +1;  % note +1 for Matlab index

%probe values of these candidates
%when probed with vector having 1 in pm_th pos
s1 = abs(prout(pm,candix));

if m>1
   
    % use vector with ones in current pos plus previous pos
    pmprev = perm(m-1);
    evpch = zeros(M,1);  evpch([pm pmprev]) = 1;
    prout_pch = probe_ye(y,evpch);
    ixpch = rec(m-1).candix-1;  % note -1 coz it was a matlab index.
    ixvec = bitxor(candix-1,ixpch);
    s2 = abs(prout_pch(ixvec+1))';
        
else
    s2 = zeros(size(s1));
    %pmprev = [];
end

% order candidates based on scores
[foo,rix] = sort(-(s1+s2));    % score + parity score
    
rmax = s1(rix);
rmax_pch = s2(rix);

% this is used in b computation.
% Set up here to avoid doing in every loop
a = (0:2^M-1)';
aval = double( dec2bin(a,M)=='1' );  % binary expansion

% process each candidate in order
Pout = [];  
crec = [];
bout = [];

for v = 1:min(nfree,params.tree_order)

    vsc = rmax(v);
    vix = rix(v);
    psc = rmax_pch(v);

    newrow = cands(vix,:);
    Phat = P;
    Phat(pm,:) = newrow;
    Phat(:,pm) = Phat(pm,:);
    
    newrec.candix = candix(vix);
    newrec.vsc = vsc;
    newrec.psc = psc;

    recx = [rec; newrec];

    %have we filled those rows of matrix which are constrained by
    %the D-G basis?
    if DGto==(DGblen)
        
        Psynth = Phat;
        
        % "detection"
        % compute Hadamard probe on the data for each row of P matrix.
        % accumulate Hadamard power and perform Sequential test until
        % decision on presence/absence of component is made.
        hadpow = zeros(2^M,1);
        detected = 0;
        Mvec = 0:2^M-1;
        for r = 1:M
            % probe data with this row of matrix.
            ev = 2^(r-1);
            evb = abs(dec2bin(ev,M)=='1');
            [hx,yp] = probe_ye(y,ev);
            % permute wrt P row, so peak should be at index zero
            prows = mod( sum(Psynth(M-r+1,:),1), 2);
            % note M-r+1 should really be find(ev)
            % but I can cheat 'coz I know location of only 1-value
            rval = (2.^(M-1:-1:0)) * prows';
            ixvec = bitxor(Mvec,rval);
            hxs = hx(ixvec+1);
            % add on this power
            hadpow = hadpow + real(hxs.*conj(hxs));
            foo = findoutliers(hadpow, 2.5);
            if foo(1)>0 % was ~=0 but we want big outliers only
                detected = 1;
            end
            
            if detected~=0
                break;  % decision made so stop loop
            end
        end 
        if detected==1
            % confirm detection and find b 
            % via HadXform 
            Phat = Psynth;
            euse = zeros(M,1); 
            Pe = rem( Phat*euse ,2);
            phi = gen_chirp(Phat,Pe');   
            foo = conj(phi).*y.*( (-1).^( aval*Pe) );
            Hfoo = fhtnat(foo);  % * scale 2^M
            Hfoo = Hfoo * 2^M;  % restore scale
            [qq,ww] = max(abs(Hfoo));
            bhat = aval(ww,:);
            
            % now test this b.
            % power in max of Had xform
            pmax = qq^2;
            % power in rest
            ptot = real(Hfoo'*Hfoo);
            prest = ptot-pmax;
            % power per bin
            prestbin = prest/(2^M-1);
            % would expect this to be power of bin sans component.
            % and thus standard deviation of amplitude should be
            sigamp = sqrt(prestbin);       
            confirmed = qq>3*sigamp;
 
            %accept candidate only if b test is good.
            if confirmed
                Pout = Phat;
                crec = recx;
                bout = bhat;
                break;  % terminate loop;  accept this candidate.
            end
        end
    else
        
        % recursion for next row
        [Pout_, bout_, crec_] = fPN_recursive(re,m+1,y,Phat,recx,params);
        
        if ~isempty(Pout_)
            % if we get here, return cand deemed OK
            Pout = Pout_;
            bout = bout_;
            crec = crec_;
            break;  % quit loop & return
        end

    end % if DGto==DGblen 
end  % for v
return;