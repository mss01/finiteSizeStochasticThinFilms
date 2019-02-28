function [h_min_rim h_centre_Joye t_rim v_re_Joye dhdt_rim dhdt_centre c_r_Joye ratio_v_vre ratio_vc_vre] = ...
                    joyeAnalysis(hJoyeStart, hJoyeEnd, h_min, h_max_dimp_l, h_max_dimp_r, h_centre_j, deltaT, seN, t_store, kappa, L_flat, R_f, h0_init, Rc);

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0, 'defaulttextInterpreter','latex');

% matFilesToRead = dir('*.mat'); 
% for q = 1:length(matFilesToRead)
%     load(matFilesToRead(q).name); 
% end 


h_min_rim = h_min(h_min < hJoyeStart & h_min > hJoyeEnd);
h_max_dimp_l_rel = h_max_dimp_l(h_min < hJoyeStart & h_min > hJoyeEnd);
h_max_dimp_r_rel = h_max_dimp_r(h_min < hJoyeStart & h_min > hJoyeEnd);
h_centre_Joye = h_centre_j(h_min < hJoyeStart & h_min > hJoyeEnd)';
t = t_store;
hfig1 = figure;
plot(t(2:end), h_min)
hold on
plot(t(2:end), h_centre_j)
legend('$h_{min}$','$h_{centre}$')
xlabel('$time$')
ylabel('$h_{min}$, $h_{centre}$')
set(gca,'FontSize',14)

set(hfig1,'Units','Inches');
pos = get(hfig1,'Position');
set(hfig1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hfig1,'h_min_h_centreVst','-dpdf','-r300')

t_rim = t(h_min < hJoyeStart & h_min > hJoyeEnd);
v_re_Joye = 16*h_min_rim.^3*kappa./L_flat^2;
v_re_Joye_centre = 16*h_centre_Joye.^3*kappa./L_flat^2;
for i = 1:length(h_min_rim)-1
    dhdt_rim(i) = (h_min_rim(i+1) - h_min_rim(i))./(seN*deltaT);
end

for i = 1:length(h_min_rim)-1
    dhdt_centre(i) = (h_centre_Joye(i+1) - h_centre_Joye(i))./(seN*deltaT);
end

c_r_Joye = 2.*h_min_rim.*h0_init.*Rc./R_f^2;
c_r_Joye_l = 2.*h_max_dimp_l_rel.*h0_init.*Rc./R_f^2;
c_r_Joye_r = 2.*h_max_dimp_r_rel.*h0_init.*Rc./R_f^2;
c_r_Joye_centre = 2.*h_centre_Joye.*h0_init.*Rc./R_f^2;
ratio_v_vre = abs(dhdt_rim)./v_re_Joye(2:end);
ratio_vc_vre = abs(dhdt_centre)./v_re_Joye_centre(2:end);
hfig2 = figure;
loglog(c_r_Joye(2:end), ratio_v_vre)
xlabel('$C_r$')
ylabel('$v/v_{re}$')
set(gca,'FontSize',14)

set(hfig2,'Units','Inches');
pos = get(hfig2,'Position');
set(hfig2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hfig2,'thinningRates_vs_C_R_hmin','-dpdf','-r300')

hfig3 = figure;
loglog(c_r_Joye_l(2:end), ratio_v_vre, 'o')
hold on
loglog(c_r_Joye_r(2:end), ratio_v_vre, 'o')
hold on
loglog(c_r_Joye_centre(2:end), ratio_vc_vre,'--')
legend('rim thinning rate-l','rim thinning rate-r','centre thinning rate','location','best')
xlabel('$C_r$')
ylabel('$v/v_{re}$')
set(gca,'FontSize',14)

set(hfig3,'Units','Inches');
pos = get(hfig3,'Position');
set(hfig3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hfig3,'thinningRates_vs_C_R01','-dpdf','-r300')


save('JoyeComparison.mat')

end