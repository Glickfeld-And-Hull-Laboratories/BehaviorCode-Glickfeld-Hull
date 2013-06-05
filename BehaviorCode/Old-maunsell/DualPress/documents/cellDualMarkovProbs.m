
% consts
pLeft = 0.1;
pSwitch = 0.0;

odds = pLeft ./ (1-pLeft);
odds = 1/odds;
b = pSwitch / (1+odds.^2) * (1+odds);
a = b * odds

fprintf(1, 'Pleft, Pswitch: %g %g;  Odds, a, b: %g %g %g\n', ...
    chop(pLeft,3), chop(pSwitch,3), ...
    chop(odds,3), chop(a,3), chop(b,3));
    

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
    trIsLeft(iS) = nextTrL == 1;
end
    
pLHat = sum(trIsLeft) ./ nTrs;
pRHat = 1-pLHat;

d0 = diff(trIsLeft);
pStayHat = sum(d0==0) ./ nTrs;
pSwitchHat = 1-pStayHat;

fprintf(1, 'est: Pleft, Pswitch: %g %g\n', ...
    chop(pLHat,3), chop(pSwitchHat,3))


%% run fn
lList = 0:0.1:1;
swList = 0:0.1:1;

nL = length(lList);
nSw = length(swList);
lMat = repmat(NaN, [nL nSw]);
swMat = repmat(NaN, [nL nSw]);

for iL = 1:nL
    for iSw = 1:nSw
        tL = lList(iL);
        tSw = swList(iL);
        
        [lHat swHat] = fnSimMarkovProb(tL, tSw, false);
        
        lMat(iL,iSw) = lHat;
        swMat(iL,iSw) = swHat;
    end
end

%% figurize
idealL = repmat(lList', [1 nSw]);
idealSw = repmat(swList, [nL 1]);

pctErrL = (lMat - idealL) ./ (idealL+0.005) * 100;
pctErrSw = (swMat - idealSw) ./ (idealSw+0.005) * 100;
6
figH = figure(3); clf; hold on;
%surf(lList, swList, lMat-idealL);
[cs h] = contour(lList, swList, swMat - idealSw);
clabel(cs,h);
xlabel('P_{left} requested');
ylabel('P_{sw} requested');
title('asymptotic Pswitch error ($\hat{P_{sw}} - P{sw}$)', 'Interpreter', 'latex');
box on; set(gca, 'TickDir', 'out');
axis square;

figH = figure(4); clf; hold on;
%surf(lList, swList, lMat-idealL);
[cs h] = contour(lList, swList, lMat - idealL);
clabel(cs,h);
xlabel('P_{left} requested');
ylabel('P_{sw} requested');
title('asymptotic Pleft error ($\hat{P_{left}} - P_{left}$)', 'Interpreter', 'latex');
box on; set(gca, 'TickDir', 'out');
axis square;

%% export
clipclip([3 4], 9*[1 0.75], 'pdf');