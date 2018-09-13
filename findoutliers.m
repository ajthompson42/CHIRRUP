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

function out = findoutliers(dat,mth)

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
