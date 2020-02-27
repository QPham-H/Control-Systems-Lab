function [Mp,tr,ts] = par(z,wn)

%Mp
if z >= 1
    Mp = 0;
else
    Mp = exp(-pi*z/sqrt(1-z*z));
end
% tr
if z < 1.2
    tr = (1.2-0.45*z+2.6*z*z)/wn;
else
    tr = (4.7*z-1.2)/wn;
end
%ts
if z > 0.69
    ts = (6.6*z-1.6)/wn;
else
    ts = -(0.5/(z*wn))*log((1-z*z)/400);
end
end