clc,clear,close all
load('bistatic_pos.mat')
load('monostatic_pos.mat')
load('mimo_pos.mat')
%%
% target position
xt = -10;
yt = 16;

D = 5;% position of the wall on the x-axis

% pos of the image
xi = -(xt-2*D);
yi = yt;

figure;
plot(0, 0, 'y-o');
grid on
hold on
axis equal
plot(xt, yt,'m*')
xline(D)
plot(xi, yi,'c^')

plot(monost_pos(1),monost_pos(2),'rx')
plot(bistat_pos(1),bistat_pos(2),'gx')
plot(MIMO_pos(1),MIMO_pos(2),'bx')

legend('array','target', 'refr. surf.', 'image', 'mono','bist','MIMO')