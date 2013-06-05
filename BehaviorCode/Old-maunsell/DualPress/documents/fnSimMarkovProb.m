function [pLeftHat, pSwitchHat] = f(pLeft, pSwitch, doPrint)

if nargin < 3, doPrint = true; end

% consts
%pLeft = 0.1;
%pSwitch = 0.0;

odds = pLeft ./ (1-pLeft);
odds = 1/odds;
b = pSwitch / (1+odds.^2) * (1+odds);
a = b * odds;

if doPrint
    fprintf(1, 'Pleft, Pswitch: %g %g;  Odds, a, b: %g %g %g\n', ...
        chop(pLeft,3), chop(pSwitch,3), ...
        chop(odds,3), chop(a,3), chop(b,3));
end

    

% sample some
nTrs = 100000;
trIsLeft = repmat(NaN, [1 nTrs]);
trIsLeft(1) = 1;
for iS = 2:nTrs
    lastTrL = trIsLeft(iS-1);
    if lastTrL == 1
        nextTrL = rand < (1-a);
    else
        nextTrL = rand < b;
    end
    
    % every 50 trs pick one unconditionally
    if mod(iS,50) == 0
        nextTrL = rand < pLeft;
    end

    trIsLeft(iS) = nextTrL == 1;
end
    
pLHat = sum(trIsLeft) ./ nTrs;
pRHat = 1-pLHat;

d0 = diff(trIsLeft);
pStayHat = sum(d0==0) ./ nTrs;
pSwitchHat = 1-pStayHat;

if doPrint
    fprintf(1, 'est: Pleft, Pswitch: %g %g\n', ...
        chop(pLHat,3), chop(pSwitchHat,3));
end

pLeftHat = pLHat;