close all

k = .0555556;
br = 198;
G = tf(k*br,[1,k*br]);
[y,t] = step(G);

ssVal = y(end);
t_10index = find(  y > .1*ssVal, 1, 'first');
% I've done the 10% index, you do the 90%:
t_90index  = find(  y > .9*ssVal, 1, 'first');
tr = t(t_90index )-t(t_10index );

plot(t,y)

