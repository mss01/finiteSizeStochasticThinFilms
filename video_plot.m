function M = video_plot(L_flat, x, Y, t_new, h0_init, t, t_scale, l_scale, beginDrainageTime_right, endDrainageTime_right, x_dimple_loc_right, res_limit, deltaX)

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0, 'defaulttextInterpreter','latex');
area(x'*l_scale*10^6,Y*h0_init*10^9);
ylim([0 1.4*h0_init*10^9]);
xlim([0 1.5*L_flat*l_scale*10^6]);
yli = yline(h0_init*10^9,'-r','h_o','LineWidth',1.5);
yli.LabelHorizontalAlignment = 'center';
pat1 = line([(x_dimple_loc_right - res_limit/2)*l_scale*10^6 (x_dimple_loc_right + res_limit/2)*l_scale*10^6],[2.5 2.5]);
pat1.Color = 'r';
pat1.LineWidth = 2;
xlabel('r ($\mu$m)','Fontsize',16);
ylabel('h (nm)','Fontsize',16);
set(gca,'FontSize',18);
if t >= beginDrainageTime_right*t_scale
    yline(100,'-.k','LineWidth',1.5);
end
if t >= endDrainageTime_right*t_scale
    yline(25,':k','LineWidth',1.5);
end
legend(t_new);
M = print('-RGBImage',sprintf('-r%d',150)); %% 150 = resolution, normally gives a good animation, and also occupies less space
end