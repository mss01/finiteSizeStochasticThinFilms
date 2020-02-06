function [h_min_rim h_centre_Joye t_rim v_re_Joye dhdt_rim dhdt_centre c_r_Joye ratio_v_vre ratio_vc_vre h_centre_FM t_series_rel h_fit_MTR bb] = ...
                    joyeAnalysis(filmConfiguration, disjPress_switch, hJoyeStart, hJoyeEnd, h_min, h_min_fit_AtDimple, h_max_dimp_r, h_centre_j, deltaT, seN, t_store, kappa, L_flat, R_f, h0_init, Rc);

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0, 'defaulttextInterpreter','latex');

% matFilesToRead = dir('*.mat'); 
% for q = 1:length(matFilesToRead)
%     load(matFilesToRead(q).name); 
% end 


%% fitting routine added
min_h_Rf_rel = h_min(30e-9/h0_init < h_min & h_min < 80e-9/h0_init);
t_series_rel = t_store(30e-9/h0_init < h_min & h_min < 80e-9/h0_init);

ff = @(bb,xx) bb(1).*exp(bb(2).*t_series_rel);                                              % Objective Function
bb = fminsearch(@(bb) norm(min_h_Rf_rel - ff(bb,t_series_rel)), [80e-9/h0_init; -0.001])                  % Estimate Parameters

% tbl = table(t_series_rel', min_h_Rf_rel');
% % modelfunc = @(bbb,xxx) bbb(1)*xxx(:,1).^bbb(2);
% modelfunc = @(bbb,xxx) bbb(1)*exp(-bbb(2).*xxx(:,1));
% beta0 = [0.5, 0.5];
% mdl = fitnlm(tbl, modelfunc, beta0);
% coeffic = mdl.Coefficients{:, 'Estimate'};
% aplu = coeffic(1);
% taplu = coeffic(2);
% h_fit_MTR = aplu*exp(-taplu.*t_series_rel);

h_fit_MTR = ff(bb,t_series_rel);

%% fitting routine ended above
h_min_rim = h_min(h_min < hJoyeStart & h_min > hJoyeEnd);
h_max_dimp_r_rel = h_max_dimp_r(h_min < hJoyeStart & h_min > hJoyeEnd);
h_centre_Joye = h_centre_j(h_min < hJoyeStart & h_min > hJoyeEnd)';
t = t_store;

h_centre_FM = (0.0004.*L_flat.^6./t).^0.25;
hfig01 = figure;
hfig01.Renderer = 'Painters';
plot(t, h_centre_FM);
xlabel('$t$');
ylabel('$h_{FM}$');
set(gca,'FontSize',14)

set(hfig01,'Units','Inches');
pos = get(hfig01,'Position');
set(hfig01,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(hfig01,'h_FMVst','-dpdf','-r300');
savefig('h_FMVst')

hfig1 = figure;
plot(t(2:end), h_min);
hold on
plot(t(2:end), h_centre_j);
legend('$h_{min}$','$h_{c}$');
xlabel('$t$');
ylabel('$h_{min}$, $h_{c}$');
set(gca,'FontSize',14);

set(hfig1,'Units','Inches');
pos = get(hfig1,'Position');
set(hfig1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(hfig1,'h_min_h_centreVst','-dpdf','-r300');
savefig('h_min_h_centreVst')

hfig1a1 = figure;
plot(t(1:end), [1 h_min]);
legend('$h_{min}$');
xlabel('$t$');
ylabel('$h_{min}$');
set(gca,'FontSize',14);

set(hfig1a1,'Units','Inches');
pos = get(hfig1a1,'Position');
set(hfig1a1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(hfig1a1,'h_min_Vst','-dpdf','-r300');
savefig('h_minVst')

hfig1a1_log = figure;
plot(t(1:end), log([1 h_min]),  'o', 'color', 'k');
hold on
plot(t_series_rel, log(h_fit_MTR), 'color', 'r', 'LineStyle', '--', 'Linewidth', 1);

% legend('$h_{min}$');
xlabel('$t$ [-]');
ylabel('$h_{min}$ [-]');
set(gca,'FontSize',14);

set(hfig1a1_log,'Units','Inches');
pos = get(hfig1a1_log,'Position');
set(hfig1a1_log,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(hfig1a1_log,'h_minLOG_Vst','-dpdf','-r300');
savefig('h_minLOGVst')

hfig1a2 = figure;
loglog(t(1:end), [1 h_min], 'o', 'color', 'k');
xlim([0 500])
ylim([5*10^-3 1])
% legend('$h_{min}$');
set(gca, 'XTick', [10^-1 1 10 100] );
xlabel('$\tilde{t}$ [-]');
ylabel('$\tilde{h}_{min}$ [-]');
set(gca,'FontSize',14);

set(hfig1a2,'Units','Inches');
pos = get(hfig1a2,'Position');
set(hfig1a2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(hfig1a2,'h_min_Vst_loglog','-dpdf','-r300');
savefig('h_minVst_loglog')

hfig1a3 = figure;
tvdW = (t_store(end) + deltaT - t_store);
loglog(tvdW(1:end), [1 h_min], '.');
legend('$h_{min}$');
xlabel('$t$');
ylabel('$h_{min}$');
set(gca,'FontSize',14);

set(hfig1a3,'Units','Inches');
pos = get(hfig1a3,'Position');
set(hfig1a3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(hfig1a3,'h_min_Vst_vdW','-dpdf','-r300');
savefig('h_minVst_vdW')

hfig1a4 = figure;
semilogy(t(1:end), [1 h_min],'color', 'k', 'Linewidth', 1);
hold on
semilogy(t_series_rel, h_fit_MTR, 'color', 'r', 'LineStyle', '--', 'Linewidth', 1);
% legend('$h_{min}$');
xlabel('$t$ [-]');
ylabel('$h_{min}$ [-]');
set(gca,'FontSize',14);

set(hfig1a4,'Units','Inches');
pos = get(hfig1a4,'Position');
set(hfig1a4,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(hfig1a4,'h_min_Vst_loglin','-dpdf','-r300');
savefig('h_minVst_loglin')



hfig1a = figure;
plot(t(2:end), h_centre_j);
hold on
plot(t, h_centre_FM);
legend('$h_{c}$', '$h_{FM}$');
xlabel('$t$');
ylabel('$h_{c}$, $h_{FM}$');
set(gca,'FontSize',14);

set(hfig1a,'Units','Inches');
pos = get(hfig1a,'Position');
set(hfig1a,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(hfig1a,'h_FM_h_centreVst','-dpdf','-r300');
savefig('h_FM_h_centreVst')

% hfig1b = figure;
% plot(t(2:end), h_max);
% hold on
% plot(t, h_centre_FM);
% legend('$h_{c}$', '$h_{FM}$');
% xlabel('$t$');
% ylabel('$h_{c}$, $h_{FM}$');
% set(gca,'FontSize',14);
% 
% set(hfig1b,'Units','Inches');
% pos = get(hfig1b,'Position');
% set(hfig1b,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
% print(hfig1b,'h_FM_h_maxVst','-dpdf','-r300');
% savefig('h_FM_h_maxVst')

t_rim = t(h_min < hJoyeStart & h_min > hJoyeEnd);

switch filmConfiguration
    case 'finiteSizedNonFlatFilms'
        %% the following is based on cartesian coordinates and accounting for both Laplace pressure and disjoining pressure
        if isequal(disjPress_switch, 'on') 
            v_re_Joye = h_min_rim.^3./L_flat^2.*(6*kappa + 1./h_min_rim.^3);
            v_re_Joye_centre = h_centre_Joye.^3./L_flat^2.*(6*kappa + 1./h_centre_Joye.^3);
        else isequal(disjPress_switch, 'off') 
            v_re_Joye = h_min_rim.^3./L_flat^2.*(6*kappa);
            v_re_Joye_centre = h_centre_Joye.^3./L_flat^2.*(6*kappa);
        end 
    case 'axisSymmetricFilm'
        %% the following is based on cylindrical coordinates and accounting for both Laplace pressure and disjoining pressure
        if isequal(disjPress_switch, 'on') 
            v_re_Joye = 16*h_min_rim.^3./L_flat^2.*(1 + 1./(6*kappa*h_min_rim.^3));                     
            v_re_Joye_centre = 16*h_centre_Joye.^3./L_flat^2.*(1 + 1./(6*kappa*h_centre_Joye.^3));
        else isequal(disjPress_switch, 'off') 
            %% this is based on cylindrical coordinates and for only Laplace pressure
            v_re_Joye = 16*h_min_rim.^3./L_flat^2;   
            v_re_Joye_centre = 16*h_centre_Joye.^3./L_flat^2;
        end     
end
%% the following is based on cylindrical coordinates and accounting for both Laplace pressure and disjoining pressure based on the earlier scaling with dominant wavelength

% v_re_Joye = h_min_rim.^3./L_flat^2.*(6*kappa + 1./h_min_rim.^3);
% v_re_Joye_centre = h_centre_Joye.^3./L_flat^2.*(6*kappa + 1./h_centre_Joye.^3);
for i = 1:length(h_min_rim)-1
    dhdt_rim(i) = (h_min_rim(i+1) - h_min_rim(i))./(seN*deltaT);
end

for i = 1:length(h_min_rim)-1
    dhdt_centre(i) = (h_centre_Joye(i+1) - h_centre_Joye(i))./(seN*deltaT);
end

c_r_Joye = 2.*h_min_rim.*h0_init.*Rc./R_f^2;
c_r_Joye_r = 2.*h_max_dimp_r_rel.*h0_init.*Rc./R_f^2;
c_r_Joye_centre = 2.*h_centre_Joye.*h0_init.*Rc./R_f^2;
ratio_v_vre = abs(dhdt_rim)./v_re_Joye(2:end);
ratio_vc_vre = abs(dhdt_centre)./v_re_Joye_centre(2:end);
hfig2 = figure;
loglog(c_r_Joye(2:end), ratio_v_vre,'o')
xlabel('$C_r$')
ylabel('$v/v_{re}$')
set(gca,'FontSize',14)
hold on
loglog(c_r_Joye_centre(2:end), ratio_vc_vre,'o')
xlabel('$C_r$')
ylabel('$v/v_{re}$')
set(gca,'FontSize',14)

set(hfig2,'Units','Inches');
pos = get(hfig2,'Position');
set(hfig2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hfig2,'thinningRates_vs_C_R_hmin','-dpdf','-r300')
savefig('thinningRates_vs_C_R_hmin')

hfig3 = figure;
loglog(c_r_Joye_r(2:end), ratio_v_vre, 'o')
hold on
loglog(c_r_Joye_centre(2:end), ratio_vc_vre,'--')
legend('rim thinning rate','centre thinning rate','location','best')
xlabel('$C_r$')
ylabel('$v/v_{re}$')
set(gca,'FontSize',14)

set(hfig3,'Units','Inches');
pos = get(hfig3,'Position');
set(hfig3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hfig3,'thinningRates_vs_C_R01','-dpdf','-r300')
savefig('thinningRates_vs_C_R01')

% for i = 1:length(h_min)-1
%     dhdt_min_all(i) = (h_min(i+1) - h_min(i))./(seN*deltaT);
% end

for i = 3:length(h_min)-2
    dhdt_min_all(i) = (h_min(i+2) - h_min(i-2))./(4*seN*deltaT);
end
size(dhdt_min_all)
size(t)
dhdt_min_all(1) = dhdt_min_all(3);
dhdt_min_all(2) = dhdt_min_all(3);
for i = 1:length(h_min)-1
    dhdt_c_all(i) = (h_centre_j(i+1) - h_centre_j(i))./(seN*deltaT);
end

hfig4 = figure;
hfig4.Renderer = 'Painters';
% size(t)
% size(dhdt_min_all)
% size(dhdt_c_all)
loglog(t(2:end-2), abs(dhdt_min_all), 'o')
xlabel('$t$')
ylabel('$|dh_{min}/dt|$')
set(gca,'FontSize',16)

set(hfig4,'Units','Inches');
pos = get(hfig4,'Position');
set(hfig4,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hfig4,'thinningRates_rim_vs_time','-dpdf','-r300')
savefig('thinningRates_rim_vs_time')

hfig5 = figure;
hfig5.Renderer = 'Painters';

loglog(t(2:end-1), abs(dhdt_c_all), 'o')
xlabel('$t$')
ylabel('$|dh_{c}/dt|$')
set(gca,'FontSize',16)

set(hfig5,'Units','Inches');
pos = get(hfig5,'Position');
set(hfig5,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hfig5,'thinningRates_c_vs_time','-dpdf','-r300')
savefig('thinningRates_c_vs_time')

hfig6 = figure;
hfig6.Renderer = 'Painters';
t_ls = (t(end) - t);
loglog(t_ls(2:end-2)./(144*kappa^2), abs(dhdt_min_all), 'o')
hold on
h_fit_ls = (0.08734/kappa^(2/5).*t_ls(2:end-1).^(-4/5));
loglog(t_ls(2:end-1)./(144*kappa^2), h_fit_ls)
xlabel('$(t_r - t)/144 \kappa^2$')
ylabel('$|dh_{min}/dt|$')
set(gca,'FontSize',16)

set(hfig6,'Units','Inches');
pos = get(hfig6,'Position');
set(hfig6,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hfig6,'thinningRates_min_vs_time_lateStage','-dpdf','-r300')
savefig('thinningRates_min_vs_time_lateStage')


close all;
save('JoyeComparison.mat')

end