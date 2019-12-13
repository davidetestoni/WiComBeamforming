function visualize2D(gc, tp,ula)


tgt_r = tp.mobileRange;
tgt_th = tp.mobileAngle(1);
int_r = tp.interfRange;
int_th = tp.interfAngle(1);


%hNull = polar(pi, 2* max_rho, 'sw');
%hPlotAxes = get(hNull, 'Parent');

%hold(hPlotAxes, 'on')

hTarget = polarplot(deg2rad(tgt_th), tgt_r, 'sg');
set(hTarget, 'MarkerFaceColor', 'r')
hold on
hInterferer = polarplot(deg2rad(int_th), int_r, 'or');
set(hInterferer, 'MarkerFaceColor', 'b')
%base station
bs = polarplot(0,0, 'dc');
set(bs, 'MarkerFaceColor', 'b')

%hold(hPlotAxes, 'off')

legend('Target','Interference','Base Station', ...
    'Location','northeastoutside');
hold off


