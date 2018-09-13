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

function [P,b] = bstr2Pb(M,bstr)

MM = ((M^2+3*M)/2); %    MMold = ((M^2+M)/2);
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

ixmid = (M+1)*M/2;  %  ixmidold = (M-1)*M/2;

% fprintf('MM %i (prev %i)   ix %i (prev %i)   len %i\n',...
%     MM,MMold,ixmid,ixmidold,length(bvec));

P = false(M);  
ix = find(tril(ones(M),0));  % indices of lower tri, including diag
P(ix) = bvec(1:ixmid);
P = P | P';  % make symmetric
b = bvec(ixmid+1:end);  

return
