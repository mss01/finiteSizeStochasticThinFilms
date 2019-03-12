function postProcessFiniteSizedFilms()

h_dimensionless = [0.05:0.05:1]';  % dimensionless film thickness vector --> will be used in the calculation of theoretical velocities

matFilesToRead = dir('*.mat'); 
for q = 1:length(matFilesToRead)
    load(matFilesToRead(q).name); 
end 

% turn everything into dimensional quantities
t_rupt_det = t_rupt.*t_scale;
t_drain_det_wholeFilm = drainageTime.*t_scale;
t_drain_det_left = drainageTime_left.*t_scale;
t_drain_det_right = drainageTime_right.*t_scale;
v_thin_min_det = abs(avg_cr_thinningRate_fit).*h0_init*10^10./t_scale;
h_cr_det_final = h_cr_final.*h0_init*10^10;

% Reynolds theory

[v_re_det t_re v_MTR t_MTR v_MTR_1997Paper v_MTR_Tsekov] = Reynolds_and_MTR(h_dimensionless, kappa, L_flat, R_f, h0_init,...
                                                                    t_scale, h_drain_start, h_drain_end, visc, gam, Rc, A_vw);

% read Radoev data
[L_film_Radoev L_film_Radoev_sort v_thin_Radoev_sort std_v_thin_Radoev h_cr_Radoev_sort std_h_cr_Radoev] = readRadoevData();
% read Manev1984 data
[R_manev t_drain_1stData std_t_drain_1stData t_drain_2ndData std_t_drain_2ndData] = Manev1984();
%% plots

% first we look at the rupture times for different film sizes
h1a = figure;
h1a.Renderer = 'Painters';
figureName_tr = strcat('ruptureTimes_','h0_',num2str(h0_init*10^9),'nm');
loglog(R_film, t_rupt_det, 'o')
xlabel('$R_{film}$ ($\mu$m)','Fontsize',14)
ylabel('$t$ (s)','Fontsize',14)
set(gca,'FontSize',14)

set(h1a,'Units','Inches');
pos = get(h1a,'Position');
set(h1a,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h1a,figureName_tr,'-dpdf','-r300')

% next we look at the drainage times for different film sizes
h1b = figure;
h1b.Renderer = 'Painters';
figureName_dr = strcat('drainageTimes_','h0_',num2str(h0_init*10^9),'nm');
errorbar(R_manev, t_drain_1stData, std_t_drain_1stData, 'd')
set(gca, 'XScale','log')
set(gca, 'YScale','log')
hold on
errorbar(R_manev, t_drain_2ndData, std_t_drain_2ndData, 'd')
set(gca, 'XScale','log')
set(gca, 'YScale','log')
hold on
loglog(R_film, t_drain_det_left, 'o')
xlabel('$R_{film}$ ($\mu$m)','Fontsize',14)
ylabel('$t_{drain}$ (s)','Fontsize',14)
hold on 
loglog(R_film, t_re, 'o')
hold on
loglog(R_film, t_MTR,'o')
set(gca,'FontSize',16)
ylim([0.1*min(t_drain_det_left) 10*max(t_drain_det_left)])
xlim([10 1100])
legend('Manev - DS1', 'Manev - DS2','$\theta$ = 0','Reynolds theory','MTR theory','Location','best')

set(h1b,'Units','Inches');
pos = get(h1b,'Position');
set(h1b,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h1b,figureName_dr,'-dpdf','-r300')

% now on to the critical thinning rates
h2b = figure;
h2b.Renderer = 'Painters';
figureName_thRates = strcat('criticalVel_','h0_',num2str(h0_init*10^9),'nm');

loglog(R_film, v_thin_min_det, 'o')
hold on
errorbar(L_film_Radoev_sort*10^6, v_thin_Radoev_sort, std_v_thin_Radoev, 'o')
set(gca, 'Xscale', 'log')
set(gca, 'Yscale', 'log')
hold on
loglog(R_film, v_re_det(:,7),'o')
hold on
loglog(R_film, v_MTR(:,7),'o')
xlabel('$R_{film}$ ($\mu$m)','Fontsize',14)
ylabel('$V$ ($\AA/s$)','Fontsize',14)
set(gca,'FontSize',16)
ylim([0.1*min(v_thin_min_det) 10*max(v_thin_min_det)])
xlim([10 1100])
legend('thinning rate, $\theta = 0$', 'Radoev thin rate', 'Reynolds thinning rate', 'MTR thinning rate')

set(h2b,'Units','Inches');
pos = get(h2b,'Position');
set(h2b,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h2b,figureName_thRates,'-dpdf','-r300')


h2c = figure;
h2c.Renderer = 'Painters';
figureName_thRates_01 = strcat('criticalVel_linLin_','h0_',num2str(h0_init*10^9),'nm');

plot(R_film, v_thin_min_det, 'o')
hold on
errorbar(L_film_Radoev*10^6, v_thin_Radoev_sort, std_v_thin_Radoev, 'o')
hold on
plot(R_film, v_re_det(:,7),'o')
hold on
plot(R_film, v_MTR(:,7),'o')
ylim([0.1*min(v_thin_min_det) 10*max(v_thin_min_det)])
xlim([10 1100])
xlabel('$R_{film}$ ($\mu$m)','Fontsize',14)
ylabel('$V$ ($\AA/s$)','Fontsize',14)
set(gca,'FontSize',14)

legend('thinning rate, $\theta$ = 0','Radoev thinning rate', 'Reynolds thinning rate', 'MTR thinning rate')

set(h2c,'Units','Inches');
pos = get(h2c,'Position');
set(h2c,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h2c,figureName_thRates_01,'-dpdf','-r300')

% critical thicknesses

h4b = figure;
h4b.Renderer = 'Painters';
figureName_crThickness = strcat('criticalThickness_','h0_',num2str(h0_init*10^9),'nm');

loglog(R_film, h_cr_det_final(3,:), 'o')
xlabel('$R_{film}$ ($\mu$m)','Fontsize',14)
ylabel('$h_{cr}$ ($\AA$)','Fontsize',14)
set(gca,'FontSize',14)
hold on
errorbar(L_film_Radoev*10^6, h_cr_Radoev_sort, std_h_cr_Radoev, 'o')
set(gca, 'Xscale', 'log')
set(gca, 'Yscale', 'log')
set(gca,'FontSize',16)

legend('$\theta$ = 0', 'Radoev expt. data', 'Vrij','MTR','Radoev','Location','best')

set(h4b,'Units','Inches');
pos = get(h4b,'Position');
set(h4b,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h4b,figureName_crThickness,'-dpdf','-r300')

end