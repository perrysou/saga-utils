function [tec] = getVerticalTEC(fArr, phiArr, rhoArr)
f1 = fArr(1);
f2 = fArr(2);
gamma = f2^2 / (f2^2 - f1^2);
delayRho = gamma*(rhoArr(:,1) - rhoArr(:,2));
delayPhi = -gamma*(phiArr(:,1) - phiArr(:,2));
end