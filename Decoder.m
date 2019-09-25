classdef Decoder
%Class responsible for decoding chirps - contains all functons pertaining to decoding
    properties
        Y       % Y{p} is a 2^m x 2^p matrix of measurements for patch p
        r       % 0, 1 or 2; 2^r patches
        parity  % parity check codes generated in the tree encoding, leave empty if only 1 patch (r=0) is used
        re      % logical: false=complex chirps, true=real chirps
        l       % length-2^r vector of parity check digits, recommend: r=0: l=0, r=1: l=[0 15], r=2: l = [0 10 10 15]
        m       % size of chirp code (2^m is length of codewords in each slot) recommend m = 6, 7 or 8
        p       % 2^p slots, require p<=(m-r)(m-r+3)/2-1 for complex
                %                    p<=(m-r)(m-r+1)/2-1 for real
        K       % number of messages
        params  % various; see below
        patches % number of patches
    end

    methods
        function self = Decoder(Y,r,l,parity,re,m,p,K,params_in)
            addpath('utils');
            self.Y=Y;
            self.r=r;
            self.parity=parity;
            self.re=re;
            self.m=m;
            self.p=p;
            self.K=K;
            self.patches=2^r;
            self.l=l;

            %default parameter values
            self.params.alpha = 0.1; %accept components if coefficient is within alpha of 1
            self.params.circuits = 3; %5; %number of circuits in peeling decoder
            self.params.sparsity_factor = 3; %factor by which number of iterations exceeds
            %expected sparsity for a given slot
            self.params.tree_order = 1; %3; %number of candidates per row in findPb tree search
            sumpropfound = 0;

            %read in parameter inputs
            if (nargin==9)
                if isfield(params_in,'alpha')
                    self.params.alpha = params_in.alpha;
                end
                if isfield(params_in,'circuits')
                    self.params.circuits = params_in.circuits;
                end
                if isfield(params_in,'sparsity_ratio')
                    self.params.sparsity_factor = params_in.sparsity_factor;
                end
                if isfield(params_in,'tree_order')
                    self.params.tree_order = params_in.tree_order;
                end
            end


        end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [output_bits, timing] = chirrup_decode(self)

        % chirrup_decode  Performs CHIRRUP decoding
        %
        % output_bits  B x K matrix of the K B-bit messages
        % timing       running time in seconds
        %
        % AJT (12/9/18)

            global outer_recov
            warning('off','MATLAB:rankDeficientMatrix');
            flag = false;
            sumpropfound=0;

            for patch = 1:self.patches

                    tic
                    outer_recov = [];
                    count = 1;
                    Yp = self.Y{patch};
                    found = [];

                    %cycle through the slots repeatedly
                    for c = 1:self.params.circuits
                        for slot = 1:2^self.p

                            %run chirp reconstruction on given slot
                            sparsity_ratio = 3;
                            recov = self.chirp_rec(Yp(:,slot),ceil(self.params.sparsity_factor*self.K/2^(self.p-1)),slot);
                            [i j] = upper_indices(self.m, self.re);

                            if (~isempty(recov(1).P))
                                %find slot pairs for each recovered component
                                for r = 1:size(recov,2)
                                    Q=recov(r).P;
                                    for k = 1:length(i)
                                        Pvec(k) = Q(i(k),j(k));
                                    end
                                    Pbvec = [Pvec recov(r).b];
                                    already = 0;
                                    if (self.re==0)
                                        trans = Pbvec(self.m*(self.m+1)/2+self.m:-1:self.m*(self.m+1)/2+self.m-self.p+1);
                                    else
                                        trans = Pbvec(self.m*(self.m-1)/2+self.m:-1:self.m*(self.m-1)/2+self.m-self.p+1);
                                    end
                                    if outofbinary(trans)==0
                                        trans(1) = 1;
                                    end

                                    if (self.re==0)
                                        if recov(r).P(1,1)==1
                                            recov(r).P(1,1) = 0;
                                            recov(r).comps(2) = slot;
                                            recov(r).comps(1) = outofbinary(mod(trans+(intobinary(slot-1,self.p))',2))+1;
                                        else
                                            recov(r).comps(1) = slot;
                                            recov(r).comps(2) = outofbinary(mod(trans+(intobinary(slot-1,self.p))',2))+1;
                                        end
                                    else
                                        if recov(r).P(1,2)==1
                                            recov(r).P(1,2) = 0;
                                            recov(r).P(2,1) = 0;
                                            recov(r).comps(2) = slot;
                                            recov(r).comps(1) = outofbinary(mod(trans+(intobinary(slot-1,self.p))',2))+1;
                                        else
                                            recov(r).comps(1) = slot;
                                            recov(r).comps(2) = outofbinary(mod(trans+(intobinary(slot-1,self.p))',2))+1;
                                        end
                                    end

                                    %accept component if its coefficient is close to 1
                                    %if (abs(recov(r).c - 2)<self.params.alpha) || (abs(recov(r).c - 1)<self.params.alpha) %alternative condition
                                    if (real(recov(r).c)>1-self.params.alpha && real(recov(r).c)<2+self.params.alpha && abs(imag(recov(r).c))<self.params.alpha)
                                    %if 1==1%<-----------------
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
                    bit_matrix = self.find_bits(outer_recov);
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
                    timings(self.patches+1) = toc;
                    output_bits=[];
                else
                    if (self.patches>1)
                        output_bits = self.tree_decoder(processed_bits);
                    else
                        output_bits = bit_matrix';
                    end
                    kmin = min(self.K,size(output_bits,2));
                    output_bits = output_bits(:,1:kmin);
                    timings(self.patches+1) = toc;
                end

                timing = sum(timings);
        end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        function recov = chirp_rec(self,y,nitlim,slot)

        % chirp_rec     Runs the chirp reconstruction algorithm, incorporating
        %               components found in other slots.
        %
        % y        measurement vector for given slot
        % nitlim   iteration limit
        % slot     slot number
        %
        % recov    multidimensional struct storing (P,b,c) for each component
        %          e.g. for component i: recov(1).P, recov(1).b, recov(1).c
        %
        % Sina Jafarpour (14/4/08) and AJT (12/9/18)

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
                    if (self.re==0)
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
                RM = Encoder.gen_chirp(recov(r).P,recov(r).b);
                allRM = [allRM RM];
            end
            cr = allRM\y;
            y = y - allRM*cr;
            ncomp = size(recov,2);
            end

            M = log2(length(y));

            while (ncomp < nitlim && norm(y)>1e-3)

                [Phat bhat] = self.findPb(y);

                %determine component from P,b
                RM = Encoder.gen_chirp(Phat,bhat);
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
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function bits = find_bits(self,res)

            bits = [];
            OK = [];
            nOK = 0;

            K = length(res);
            s = 1;

            for k = 1:K

                Pk = res(k).P;
                bk = res(k).b;
                ck = res(k).c;
                compsk = res(k).comps;
                [i j] = upper_indices(self.m, self.re);
                for l = 1:length(i)
                    Pvec(l) = Pk(i(l),j(l));
                end
                bits2 = [Pvec bk];
                bits2 = bits2(1,2:end);
                slot_bit = (intobinary(compsk(1)-1,self.p))';
                bits(k,:) = [bits2 slot_bit];
            end
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function found_bits = tree_decoder(self, output_bits)

        % tree_decoder  Decodes a patched message encoded by tree_encoder.m
        %
        % output_bits     output_bits{i} is a B x k_i matrix of the k_i B-bit
        %                 messages found in patch i
        %
        % found_bits       N x K matrix of the K repatched N-bit messages
        %
        % Code is based on 'A Coupled Compressive Sensing Scheme for Unsourced
        % Multiple Access' by Amalladinne et al. 2018 (arXiv 1806.00138)

            n = length(output_bits);
            for i = 1:n
                k(i) = size(output_bits{i},2);
            end
            N = size(output_bits{1},1);

            found_bits = [];

            for i = 1:k(1)

                cand_bits = output_bits{1}(:,i);

                for j = 2:n

                    cand_bits = [repmat(cand_bits,[1 k(j)]); kron(output_bits{j}(:,:),ones(1,size(cand_bits,2)))];
                    parity_check = mod(self.parity{j}*cand_bits(1:end-self.l(j),:)+cand_bits(end-self.l(j)+1:end,:),2);
                    find_correct = sum(parity_check)<0.5;
                    cand_bits = cand_bits(1:end-self.l(j),find_correct);

                end

                found_bits = [found_bits cand_bits];

            end
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [Pout bout] = findPb(self,y)

        % findPb  Find P matrix & b vector recursively.
        %
        % y         signal, length 2^M
        %
        % Pout      P matrix
        % bout      b vector
        %
        % Sina Jafarpour (29/4/08) and AJT (12/9/18)

            global perm prout DGbasis;
            % so these are visible from within recursive function
            twoM = length(y);
            M = log2(twoM);

            DGbasis = makeDGC(self.re,M);

            if (self.re==0)
                DGblen = M*(M+1)/2;
            else
                DGblen = M*(M-1)/2;
            end

            %determine probe response to each basis vector
            %probe = Hadamard xform of y^*(a)y(a+e), e = basis vector
            evp = eye(M);
            prout = zeros(M,twoM);
            for m = 1:M
                prout(m,:) = self.probe_ye(y,evp(m,:));
            end

            perm = (1:M)';

            %manipulate permutation so only rows encompassing the first
            %DGblen elements will be used.
            if (self.re==0)
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
            [Pout,bout,rec] = self.fPN_recursive(1,y,zeros(M),[]);

            %what to do if nothing found?
            if isempty(Pout)
                % Oh heck.  Didn't find anything.
                % Make one up at random...  hopefully this will kick the signal enough
                % to make something stand out next time, and this bogus component will
                % be subsequently removed.
                [Pout,bout] = self.bstr2Pb(M,[]);
            end
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [Pout,bout,crec] = fPN_recursive(self,m,y,P,rec)

            global perm prout DGbasis;
            %defined prior to call of this recursive function

            M = length(P);

            %which basis matrices go with current row
            DGblen = size(DGbasis,3);
            if (self.re==0)
                csM1 = cumsum(M:-1:1);
            else
                csM1 = cumsum(M-1:-1:1);
            end
            DGto = csM1(m);
            DGto = min(DGto,DGblen);

            pm = perm(m);
            row = P(pm,:);
            if (self.re==0)
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
                prout_pch = self.probe_ye(y,evpch);
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

            for v = 1:min(nfree,self.params.tree_order)

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
                        [hx,yp] = self.probe_ye(y,ev);
                        % permute wrt P row, so peak should be at index zero
                        prows = mod( sum(Psynth(M-r+1,:),1), 2);
                        % note M-r+1 should really be find(ev)
                        % but I can cheat 'coz I know location of only 1-value
                        rval = (2.^(M-1:-1:0)) * prows';
                        ixvec = bitxor(Mvec,rval);
                        hxs = hx(ixvec+1);
                        % add on this power
                        hadpow = hadpow + real(hxs.*conj(hxs));
                        foo = self.findoutliers(hadpow, 2.5);
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
                        phi = Encoder.gen_chirp(Phat,Pe');
                        foo = conj(phi).*y.*( (-1).^( aval*Pe) );
                        Hfoo = self.fhtnat(foo);  % * scale 2^M
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
                    [Pout_, bout_, crec_] = self.fPN_recursive(m+1,y,Phat,recx);

                    if ~isempty(Pout_)
                        % if we get here, return cand deemed OK
                        Pout = Pout_;
                        bout = bout_;
                        crec = crec_;
                        break;  % quit loop & return
                    end

                end % if DGto==DGblen
            end  % for v

        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [prb,yprod] = probe_ye(self,y,ein)

        % probe_ye  probe data vector with error vector
        %
        % given data vector y and "error value" e,
        % compute Hadamard transform prb  of yprod = conj(y(a)) * y(a+e)


        % Sina Jafarpour 28/3/08

            Mpow = length(y);
            M = log2(Mpow);  % this should be an integer
            if isscalar(ein)
                e = ein;
            else
                e = (2.^(M-1:-1:0))*ein(:);
            end

            % index vector (starting at 0)
            a = (0:Mpow-1)';

            % a plus e values
            ape = bitxor(a,e);

            % y(a+e)
            yape = y(ape+1);

            % product
            yprod = conj(y).*yape;

            % use fast Walsh-Hadamard transform routine
            prb = self.fhtnat(yprod);

        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function x=fhtnat(self, data)
        %------------------------------------------------------
        %1D Natural(Hadamard)ordered Fast Hadamard Transform
        %------------------------------------------------------
        % Author: Gylson Thomas
        % e-mail: gylson_thomas@yahoo.com
        % Asst. Professor, Electrical and Electronics Engineering Dept.
        % MES College of Engineering Kuttippuram,
        % Kerala, India, February 2005.
        % copyright 2007.
        % This function implements the 1D natural(Hadamard)ordered Fast Hadamard Transform,
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
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [P,b] = bstr2Pb(self,M,bstr)
        % convert a bitstring (as an integer or as a logical vector) to P and b
        % * updated for Kerdock codes:  no longer zero diagonal.
        %
        % [P,b] = bstr2Pb(M,[]);                randomly assigned P and b
        % [P,b] = bstr2Pb(M,[1 0 1 0 ...])      take bits from vector elements
        % [P,b] = bstr2Pb(M,int);               take bits from int

        % began life as a private function of genRM.m
        % SJS 6/11/07
        % SJS 28/3/08 rehashed for Kerdock codes:
        %             diagonal no longer constrained to be all-zeros
        %             total # free elements thus M(M+3)/2

            if self.re == 0
                MM = ((M^2+3*M)/2); %    MMold = ((M^2+M)/2);
                ixmid = (M+1)*M/2;
            else
                MM = 0.5*M*(M+1);
                ixmid = (M-1)*M/2;
            end

            if isempty(bstr)
                %bstr = floor(rand*2^MM);  % precision errors!
                bstr = randn(1,MM)>0;
            end
            if isscalar(bstr)
                % could use dec2bin...
                bvec = false(1,MM);
                for q = 1:MM
                    bit = rem(bstr,2); % bitand(bstr,1);
                    bvec(end-q+1) = bit;
                    bstr = floor(bstr/2);
                end
            else
                bvec = bstr;
                % should check that length is M, but what to do if it ain't?
            end


            % fprintf('MM %i (prev %i)   ix %i (prev %i)   len %i\n',...
            %     MM,MMold,ixmid,ixmidold,length(bvec));
            P = false(M);
            if self.re == 0
                ix = find(tril(ones(M),0));  % indices of lower tri, including diag
            else
                ix = find(tril(ones(M),-1)); % doesnt include the diagonal
            end
            P(ix) = bvec(1:ixmid);
            P = P | P';  % make symmetric
            b = bvec(ixmid+1:end);

            end
        


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function out = findoutliers(self,dat,mth)

        % findoutliers  Find outliers in a single data vector.
        %
        % foo = findoutliers(dat,mth)
        %       dat must be a vector.
        %
        % foo = nonzero where dat element was an outlier
        %       where possible, sign of foo indicates "direction" of outlier
        %
        % mth = method:
        %       1   = 3 sigmas from the mean
        %       2   = quartile test, Mendenhall & Sincich quartiles
        %       2.5 = quartile test, Tukey quartiles
        %
        % see  http://mathforum.org/library/drmath/view/52720.html

            if nargin<2 | isempty(mth)
                mth = 2;
            end
            ndat = length(dat);

            switch floor(mth)
                case 1
                    % 3 standard deviations from mean
                    if max(size(dat))<=1
                        out = zeros(size(dat));
                        return
                    else
                        sd = max(std(dat),eps);
                    end

                    mu = mean(dat);
                    demu = abs(dat-mu)/sd;
                    out = demu>3;

                    % introduce sign in order to tell if greater or lesser.
                    out = out.*sign(dat-mu);

                case 2
                    [sdat,ix] = sort(dat);
                    [foo,undoix] = sort(ix);  % for unsorting
                    %pctile = (1:ndat)/ndat;

                    % compute quartiles
                    switch mth-floor(mth)
                        case 0
                            % Mendenhall and Sincich method
                            L = ceil(0.25*(ndat+1));
                            U = floor(0.75*(ndat+1));
                            LQ = sdat(L);
                            UQ = sdat(U);
                        case .5
                            % Tukey method
                            LQ = median(sdat(1:ceil(ndat/2)));
                            UQ = median(sdat(ceil((ndat+1)/2):ndat));
                    end

                    IQR = UQ-LQ;
                    out = (sdat>UQ+1.5*IQR) - (sdat<LQ-1.5*IQR);
                    out = out(undoix);
                    %fprintf('findoutliers:  LQ %g  UQ %g  IQR %g  no %g\n',LQ,UQ,IQR,sum(out~=0))

                otherwise
                    error(['no such method:  ' int2str(mth)])
            end
        end



    end
end