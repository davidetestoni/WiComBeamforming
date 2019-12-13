function visualize2D(gc, tp,ula)


%hNull = polar(pi, 2* max_rho, 'sw');
%hPlotAxes = get(hNull, 'Parent');

%hold(hPlotAxes, 'on')

htgt1 = polarplot(deg2rad(tp.tgt1_angle), tp.tgt1_range, 'sg');
set(htgt1, 'MarkerFaceColor', 'r')
hold on

htgt2 = polarplot(deg2rad(tp.tgt2_angle), tp.tgt2_range, 'og');
set(htgt2, 'MarkerFaceColor', 'r')


hintf1 = polarplot(deg2rad(tp.intf1_angle), tp.intf1_range, 'sr');
set(hintf1, 'MarkerFaceColor', 'b')

hintf2 = polarplot(deg2rad(tp.intf2_angle), tp.intf2_range, 'or');
set(hintf2, 'MarkerFaceColor', 'b')
%base station
bs = polarplot(0,0, 'dc');
set(bs, 'MarkerFaceColor', 'b')

%hold(hPlotAxes, 'off')

legend('Target','Interference','Base Station', ...
    'Location','northeastoutside');
hold off


