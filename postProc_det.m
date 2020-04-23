function [t_rupt drainageTime drainageTime_right drainageTime_right_rupt  avg_cr_thinningRate_fit h_cr_final ...
                            h_cr_final_FullFilmavg R_film_final Press vol_end area_end ratio_vol ratio_area t_series_rel h_fit_MTR bb h_c_tr] = postProc_det(filmConfiguration, disjPress_switch, correctionLP_switch, R_f, L_flat, L_curv, transitionLength, ...
                                                                    deltaX, deltaT, kappa, t_scale, l_scale, seN, animationSkip, h_drain_start, h_drain_end, ...
                                                                    h_critical_start, h_critical_end, t_cr, res_limit, hJoyeStart, hJoyeEnd, h0_init, Rc, gam);

tic

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0, 'defaulttextInterpreter','latex');

folders = dir('realization*');
load('./realization1/hData.mat');
s = 1;
start = 0;
last = 0;

%%

load('realization1/hData.mat');
t = t_store;   

%% simulation set up

switch filmConfiguration
    case 'finiteSizedNonFlatFilms'
        L = 2*(L_curv + L_flat);                    % total length of the film (curved+flat)
    case 'axisSymmetricFilm'
        L = L_curv + L_flat;                        % total length of the film (curved+flat)
end
[h x] = initialProfile(kappa,L_flat,L_curv, R_f, Rc,transitionLength,deltaX, filmConfiguration, correctionLP_switch);       % get the initial conditions
tt = seN*deltaT;                                            % the rate at which files were saved
cr_thickness = 0.627*kappa^(-2/7);                          % theoretical prediction of critical thickness  
h_store_new = [h h_store];
t_store_new = [0 t_store];


%%

hfig01 = figure;
hfig01.Renderer = 'Painters';

idxTimeStamp = findTimeSeriesBy2(t_store);

for i = 1:length(idxTimeStamp)
    dummy_01(:,i) = h_store_new(:, idxTimeStamp(i));
    dummy(:,i) = vertcat(dummy_01(1,i), dummy_01(:,i));
    x_dummy = [0; x];
    plot(x_dummy*l_scale*10^6, dummy(:,i)*h0_init*10^9,'Linewidth',1, 'color', 'k', 'LineStyle', '-'); % [0.929 0.694 0.125] yellow; [0.85 0.325 0.098] red; 
%   [0 0.447 0.741] blue
%     ylim([0 1.5*h0_init*10^9]);
    ylim([0 3000]);
%     ylim([0 1.5*h0_init*10^9]);
%     xlim([0 L_flat*l_scale*10^6 + 40]);
    xlim([0 650]);
%     xlim([0 L_flat*l_scale*10^6+100]);
%     daspect([0.1 1 1])
    hold on
end
% breakxaxis([50 950]);  % for h0 = 2000 nm and Rfilm = 4000 mu m
% breakxaxis([20 3550]);  % for h0 = 300 nm and Rfilm = 4000 mu m

xlabel('$r$ ($\mu$m)','Fontsize',16)
ylabel('$h$ (nm)','Fontsize',16)
set(gca,'XTick',[0 200 400 600]);
% set(gca,'YTick',[0 100 200 300]);
set(gca,'FontSize',16)

set(hfig01,'Units','Inches');
pos = get(hfig01,'Position');
set(hfig01,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hfig01,'heightProfiles_successiveDivideBy2','-dpdf','-r300')
savefig('heightProfiles_successiveDivideBy2')

ylim([0 4*h0_init*10^9]);
xlim([0 L_flat*l_scale*10^6 + 100]);
% ylim([0 4000]);
% xlim([0 200]);

xlabel('$r$ ($\mu$m)','Fontsize',14)
ylabel('$h$ (nm)','Fontsize',14)
% set(gca,'XTick',[0 50 100 150 200]);
set(gca,'FontSize',14)

set(hfig01,'Units','Inches');
pos = get(hfig01,'Position');
set(hfig01,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hfig01,'heightProfiles_successiveDivideBy2_yRangeLarge','-dpdf','-r300')
savefig('heightProfiles_successiveDivideBy2_yRangeLarge')


%% calculate pressure along the film


[Press CapillaryPres DisjoiningPres] = calcPress(h, h_store, h_store_new, t_store, x, L_flat, deltaX, kappa);


%% identification of the extent of the flat films, usually where the dimple forms

[locDimple_right x_centre] = keyLocationsInFilm(filmConfiguration, x, L_flat, deltaX);
% [locDimple_left locDimple_right x_centre] = keyLocationsInFilm(x, L_flat, deltaX);

%% from here on starts the averaging over multiple realizations
[h_right_avg_j h_min h_centre_j v_thin_rim v_thin_centre avg_cr_thinningRate_fit h_cr_final ...
    h_cr_final_FullFilmavg drainageTime drainageTime_right drainageTime_right_rupt t_rupt h_max_dimp_r ...
    beginDrainageTime_right endDrainageTime_right x_dimple_loc_right vol_end area_end h_min_fit_AtDimple] = ...
                                    loopingOverRealizations(x, h0_init, L_flat, transitionLength, tt, x_centre, locDimple_right, res_limit, cr_thickness, ...
                                                                    h_drain_start, h_drain_end, t_cr, deltaX, deltaT, h_critical_start, h_critical_end);
 
h_c_tr = h_centre_j(end)*h0_init*10^9;  
%%                                                              

vol_1 = 2*pi*trapz(x(1:round(L_flat/deltaX)), h(1:round(L_flat/deltaX)).*x(1:round(L_flat/deltaX)));
area_film_1 = trapz(x(1:round(L_flat/deltaX)), h(1:round(L_flat/deltaX)));

ratio_vol = vol_end/vol_1;
ratio_area = area_end/area_film_1;

R_film_t = x_dimple_loc_right*l_scale;                                                              
R_film_final = x_dimple_loc_right(end)*l_scale;
% [h_right_avg_j h_left_avg_j h_avg h_min h_centre_j v_thin_rim v_thin_centre avg_cr_thinningRate_fit h_cr_final ...
%     h_cr_final_FullFilmavg drainageTime drainageTime_left drainageTime_right drainageTime_right_rupt drainageTime_left_rupt t_rupt h_max_dimp_l h_max_dimp_r] = ...
%                                     loopingOverRealizations(x, L_flat, transitionLength, tt, x_centre, locDimple_left, locDimple_right, res_limit, cr_thickness, ...
%                                                                     h_drain_start, h_drain_end, t_cr, deltaX, deltaT, h_critical_start, h_critical_end);
    
[h_min_rim h_centre_Joye t_rim v_re_Joye dhdt_rim dhdt_centre c_r_Joye ratio_v_vre ratio_vc_vre h_centre_FM t_series_rel h_fit_MTR bb] = ...
    joyeAnalysis(filmConfiguration, disjPress_switch, hJoyeStart, hJoyeEnd, h_min, h_min_fit_AtDimple, h_max_dimp_r, h_centre_j, deltaT, seN, t_store, kappa, L_flat, R_f, h0_init, Rc);

hfig2 = figure;
hfig2.Renderer = 'Painters';

t = [0 t];
x_FM = [0:deltaX:L_flat];
% y_FM = zeros(size(h_centre_FM
y_FM = h_centre_FM'.*(1 - x_FM.^2/L_flat^2);
x_dimple_loc_right = [L_flat; x_dimple_loc_right(:,1)];

plot(t*t_scale, x_dimple_loc_right(:,1)*l_scale*10^6,'o')
xlabel('t (s)','Fontsize',14)
ylabel('$R_{film}~(\mu$m)','Fontsize',14)
set(gca,'FontSize',14)

set(hfig2,'Units','Inches');
pos = get(hfig2,'Position');
set(hfig2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hfig2,'R_film_Vs_T','-dpdf','-r300')
savefig('R_film_Vs_T')

hfig2a = figure;
hfig2a.Renderer = 'Painters';

h_max_dimp_r = [h_max_dimp_r(1) h_max_dimp_r h_max_dimp_r(end)];
plot(t*t_scale, h_max_dimp_r*h0_init*10^9, 'o')

hold on
plot(t(1:end-1)*t_scale, [1; h_centre_j]*h0_init*10^9, 'o')
xlabel('t (s)','Fontsize',14)
ylabel('$h_{max}$, $h_c$ (nm)','Fontsize',14)
legend('$h_{max}$','$h_{c}$');
set(gca,'FontSize',16)


set(hfig2a,'Units','Inches');
pos = get(hfig2a,'Position');
set(hfig2a,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hfig2a,'h_max_h_c_Vs_T','-dpdf','-r300')
savefig('h_max_h_c_Vs_T')

hfig1b = figure;
plot(t, h_max_dimp_r);
hold on
plot(t(2:end), h_centre_FM);
legend('$h_{max}$', '$h_{FM}$');
xlabel('$t$');
ylabel('$h_{max}$, $h_{FM}$');
set(gca,'FontSize',14);

set(hfig1b,'Units','Inches');
pos = get(hfig1b,'Position');
set(hfig1b,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(hfig1b,'h_FM_h_maxVst','-dpdf','-r300');
savefig('h_FM_h_maxVst')
% 
% hfig_FM = figure;
% hfig_FM.Renderer = 'Painters';
% 
% for i = 1:length(idxTimeStamp)
%     dummy_01(:,i) = h_store_new(:, idxTimeStamp(i));
%     dummy(:,i) = vertcat(dummy_01(1,i), dummy_01(:,i));
%     x_dummy = [0; x];
%     plot(x_dummy*l_scale*10^6, dummy(:,i)*h0_init*10^9,'Linewidth',1, 'color', [0 0.447 0.741], 'LineStyle', '-'); % [0.929 0.694 0.125] yellow; [0.85 0.325 0.098] red; 
%     hold on
%     plot(x_FM*l_scale*10^6, y_FM(i,:)*h0_init*10^9, 'Linewidth', 1, 'color', [0.929 0.694 0.125], 'LineStyle', '-.');
% %   [0 0.447 0.741] blue
% %     ylim([0 1.5*h0_init*10^9]);
%     ylim([0 500]);
% %     xlim([0 L_flat*l_scale*10^6 + 40]);
%     xlim([0 80]);
%    
%     hold on
% end
% % breakxaxis([50 1200]);
% xlabel('$r$ ($\mu$m)','Fontsize',14)
% ylabel('$h$ (nm)','Fontsize',14)
% set(gca,'XTick',[0 20 40 60 80]);
% set(gca,'YTick',[0 100 200 300 400 500]);
% set(gca,'FontSize',14)
% 
% set(hfig_FM,'Units','Inches');
% pos = get(hfig_FM,'Position');
% set(hfig_FM,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(hfig_FM,'heightProfiles_FMCompare','-dpdf','-r300')
% savefig('heightProfiles_FMCompare')

close all;
% delete hfig1b hfig2a hfig2 hfig1_ruptm1 hfig1_rupt hfig1_dim hfig1 hfig0 hfig01 hfig02 hfig_c

save('workspace_deterministic_t_cr.mat')
save('lastFrame01.mat','x','h_store','t_store','l_scale','t_scale','h0_init', 'R_film_final', 'R_f', 'h_centre_j')

% makeAnimation_det(filmConfiguration, correctionLP_switch, animationSkip,kappa, L_flat, L_curv, R_f, Rc, transitionLength,deltaX, h_store, t_store, h0_init, t_scale, l_scale, beginDrainageTime_right, endDrainageTime_right, x_dimple_loc_right, res_limit, h_drain_start, h_drain_end);
% makeAnimation_det(filmConfiguration, correctionLP_switch, animationSkip,kappa, L_flat, L_curv, R_f, Rc, transitionLength,deltaX, h_store, t_store, h0_init, t_scale, l_scale, beginDrainageTime_right, endDrainageTime_right, x_dimple_loc_right, res_limit, h_drain_start, h_drain_end);

toc

close all

end
