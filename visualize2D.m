function visualize2D(gc, tp,ula)


%hNull = polar(pi, 2* max_rho, 'sw');
%hPlotAxes = get(hNull, 'Parent');

%hold(hPlotAxes, 'on')

htgt1 = polarplot(deg2rad(tp.tgt1_angle), tp.tgt1_range, 'sy');
set(htgt1, 'MarkerFaceColor', 'g')
hold on

htgt2 = polarplot(deg2rad(tp.tgt2_angle), tp.tgt2_range, 'ob');
set(htgt2, 'MarkerFaceColor', 'g')


hintf1 = polarplot(deg2rad(tp.intf1_angle), tp.intf1_range, 'sy');
set(hintf1, 'MarkerFaceColor', 'r')

hintf2 = polarplot(deg2rad(tp.intf2_angle), tp.intf2_range, 'ob');
set(hintf2, 'MarkerFaceColor', 'r')
%base station
bs = polarplot(0,0, 'dc');
set(bs, 'MarkerFaceColor', 'b')

%hold(hPlotAxes, 'off')

legend('Target1','Target2','Interference1','Interference2','Base Station', ...
    'Location','northeastoutside');
hold off


