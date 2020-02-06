function M = video_plot(L_flat, x, Y, t_new, h0_init, t, t_scale, l_scale, beginDrainageTime_right, endDrainageTime_right, x_dimple_loc_right, res_limit, deltaX, h_drain_start, h_drain_end, x_FM, y_FM_rel)

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0, 'defaulttextInterpreter','latex');
x = [0; x];
Y = [Y(1); Y];

area(x'*l_scale*10^6,Y*h0_init*10^9);
% hold on
% plot(x_FM*l_scale*10^6, y_FM_rel*h0_init*10^9);
ylim([0 1.4*h0_init*10^9]);
xlim([0 ((L_flat)*l_scale*10^6 + 100)]);
yli = yline(h0_init*10^9,'--r','$h_o$','LineWidth',1.5,'Fontsize',16, 'Interpreter', 'latex');
% yli.LabelHorizontalAlignment = 'center';
% pat1 = line([(x_dimple_loc_right - res_limit(1)/2)*l_scale*10^6 (x_dimple_loc_right + res_limit(1)/2)*l_scale*10^6],[2.5 2.5]);
% pat1.Color = 'r';
% pat1.LineWidth = 2;
xlabel('r ($\mu$m)','Fontsize',16);
ylabel('h (nm)','Fontsize',16);
set(gca,'FontSize',18);


legend(t_new);

%% dimensionless videos

% area(x',Y);
% ylim([0 1.4]);
% xlim([0 1.5*L_flat]);
% yli = yline(1,'-r','h_o','LineWidth',1.5);
% yli.LabelHorizontalAlignment = 'center';
% pat1 = line([(x_dimple_loc_right - res_limit(1)/2) (x_dimple_loc_right + res_limit(1)/2)],[2.5 2.5]);
% pat1.Color = 'r';
% pat1.LineWidth = 2;
% xlabel('r [-]','Fontsize',16);
% ylabel('h [-]','Fontsize',16);
% set(gca,'FontSize',18);
% % endDrainageTime_right
% % size(endDrainageTime_right)
% if t >= beginDrainageTime_right
%     yline(h_drain_start,'-.k','LineWidth',1.5);
% end
% if t >= endDrainageTime_right(1)
%     yline(h_drain_end(1),':k','LineWidth',1.5);
% end
% if t >= endDrainageTime_right(2)
%     yline(h_drain_end(2),':k','LineWidth',1.5);
% end
% if t >= endDrainageTime_right(3)
%     yline(h_drain_end(3),':k','LineWidth',1.5);
% end
% if t >= endDrainageTime_right(4)
%     yline(h_drain_end(4),':k','LineWidth',1.5);
% end
% legend(t_new);


M = print('-RGBImage',sprintf('-r%d',150)); %% 150 = resolution, normally gives a good animation, and also occupies less space
end