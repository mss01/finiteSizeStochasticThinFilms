clear all
close all
clc
tic
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaulttextInterpreter','latex');
% set(0, 'defaultLegendInterpreter','latex');

%% last frames

% lets first plot all the last frame on top of each other


folders = dir('Rf*');
a2 = cellfun(@num2str, struct2cell(folders), 'UniformOutput', false);
Out2 = sortrows(a2.',6);

hfig1 = figure;
hfig1.Renderer = 'Painters';
figs2keep = [1];
ColOrd = get(gca, 'colororder');
ColOrd([8:19],:) = parula(12);  
for j = 1:max(size(Out2)) % max(size(Out2))
    newFolder = Out2{j};
    cd(newFolder)
    load('lastFrame01.mat')
    R_film_imposed(j)  = R_f*10^6;
    R_film_final_tr(j) = R_film_final*10^6;
    h_centre_at_tr(j)  = h_centre_j(end)*h0_init*10^9;
    Pres_h_max_01(j)   = Pres_h_max; 
%     all_figs = findobj(0, 'type', 'figure');
%     delete(setdiff(all_figs, figs2keep));
    plot(x*l_scale*10^6, h_store(:,end)*h0_init*10^9, 'color', ColOrd(j,:), 'Linewidth', 1.5);
    Legend_r{j} = strcat('$R_f = $', num2str(R_f*10^6),' $\mu$m');
    ylim([0 1500])
    hold on
    cd ..
end

h_legend_r = legend(Legend_r, 'interpreter','latex','location','best'); 

xlabel('$r$ ($\mu$m)','Fontsize',14);
ylabel('$h$ (nm)','Fontsize',14);
set(gca,'FontSize',14);

set(gcf, 'unit', 'inches');
figure_size =  get(gcf, 'position');

set(h_legend_r, 'location', 'northeastoutside');
set(h_legend_r, 'unit', 'inches');
legend_size = get(h_legend_r, 'position');
figure_size(3) = figure_size(3) + legend_size(3);
set(gcf, 'position', figure_size);

set(hfig1,'Units','Inches');
pos = get(hfig1,'Position');
set(hfig1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(hfig1,'lastRuptureEvent_diffFilmThickness_test','-dpdf','-r300');

%%

hfig7a = figure;
hfig7a.Renderer = 'Painters';

for j = 1:max(size(Out2)) % max(size(Out2))
    newFolder = Out2{j};
    cd(newFolder)
    load('hProbeAvg.mat');
    Pres_h_max_01(j)   = Pres_h_max(end);
    cd ..
end

loglog(R_film_imposed, Pres_h_max_01, 'o');

xlabel('$R_{film}$ ($\mu$m)','Fontsize',14);
ylabel('$P_{hmax}$','Fontsize',14);
set(gca,'FontSize',14);

set(hfig7a,'Units','Inches');
pos = get(hfig7a,'Position');
set(hfig7a,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(hfig7a,'Pres_dimpleVs_Rfilm','-dpdf','-r300');

%%

hfig5 = figure;
hfig5.Renderer = 'Painters';

loglog(R_film_imposed, h_centre_at_tr, 'o');

xlabel('$R_{film}$ ($\mu$m)','Fontsize',14);
ylabel('$h(r=0,t=t_r)$ (nm)','Fontsize',14);
set(gca,'FontSize',14);

set(hfig5,'Units','Inches');
pos = get(hfig5,'Position');
set(hfig5,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(hfig5,'h_centreVsR_film_loglog','-dpdf','-r300');

%%

hfig6a = figure;
hfig6a.Renderer = 'Painters';

plot(R_film_final_tr, h_centre_at_tr, 'o');

xlabel('$R_{film}$ ($\mu$m)','Fontsize',14);
ylabel('$h(r=0,t=t_r)$ (nm)','Fontsize',14);
set(gca,'FontSize',14);

set(hfig6a,'Units','Inches');
pos = get(hfig6a,'Position');
set(hfig6a,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(hfig6a,'h_centreVsR_film_at_tr','-dpdf','-r300');

%%

hfig2 = figure;
hfig2.Renderer = 'Painters';


for j = 1:max(size(Out2))  % max(size(Out2))
    newFolder = Out2{j};
    cd(newFolder)
    load('lastFrame.mat')
    plot(x*l_scale/R_film_final, h_store(:,end)*h0_init*10^9, 'color', ColOrd(j,:), 'Linewidth', 1.5);
    Legend_r2{j} = strcat('$R_f = $', num2str(R_f*10^6),' $\mu$m');
    ylim([0 1500]);
    xlim([0 2]);
    hold on
    cd ..
end
h_legend_r2 = legend(Legend_r2, 'interpreter','latex','location','best'); 


xlabel('$r$ [-]','Fontsize',14);
ylabel('$h$ (nm)','Fontsize',14);
set(gca,'FontSize',14);

set(gcf, 'unit', 'inches');
figure_size_2 =  get(gcf, 'position');

set(h_legend_r2, 'location', 'northeastoutside');
set(h_legend_r2, 'unit', 'inches');
legend_size = get(h_legend_r2, 'position');
figure_size_2(3) = figure_size_2(3) + legend_size(3);
set(gcf, 'position', figure_size_2);

set(hfig2,'Units','Inches');
pos = get(hfig2,'Position');
set(hfig2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(hfig2,'lastRuptureEvent_diffFilmThickness_scaled_R_film','-dpdf','-r300');

%% 

% hfig2a = figure;
% hfig2a.Renderer = 'Painters';
% 
% for j = 1:max(size(Out2))  % max(size(Out2))
%     newFolder = Out2{j};
%     cd(newFolder)
%     load('lastFrame.mat')
%     plot(x*l_scale/R_film_final, h_store(:,end)*h0_init*10^9, 'color', ColOrd(j,:), 'Linewidth', 1.5)
%     Legend_r3{j} = strcat('$R_f = $', num2str(R_f*10^6),' $\mu$m');
%     ylim([0 1500])
%     xlim([0 1])
%     hold on
%     cd ..
% end
% h_legend_r3 = legend(Legend_r3, 'interpreter','latex','location','best'); 
% 
% 
% xlabel('$r$ [-]','Fontsize',14)
% ylabel('$h$ (nm)','Fontsize',14)
% set(gca,'FontSize',14)
% 
% set(gcf, 'unit', 'inches');
% figure_size_3 =  get(gcf, 'position');
% 
% set(h_legend_r3, 'location', 'northeastoutside');
% set(h_legend_r3, 'unit', 'inches');
% legend_size = get(h_legend_r3, 'position');
% figure_size_3(3) = figure_size_3(3) + legend_size(3);
% set(gcf, 'position', figure_size_3);
% 
% 
% set(hfig2a,'Units','Inches');
% pos = get(hfig2a,'Position');
% set(hfig2a,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
% print(hfig2a,'lastRuptureEvent_diffFilmThickness_scaled_R_film_0_1','-dpdf','-r300');



%%

% hfig2b = figure;
% hfig2b.Renderer = 'Painters';
% 
% for j = 1:max(size(Out2))  % max(size(Out2))
%     newFolder = Out2{j};
%     cd(newFolder);
%     load('lastFrame.mat');
% %     all_figs = findobj(0, 'type', 'figure');
% %     delete(setdiff(all_figs, figs2keep));
%     plot(x*l_scale/R_film_final, h_store(:,end)*h0_init*10^9/(R_f*10^6), 'color', ColOrd(j,:), 'Linewidth', 1.5);
%     xlim([0 1])
%     Legend_r4{j} = strcat('$R_f = $', num2str(R_f*10^6),' $\mu$m');
%     hold on
%     cd ..
% end
% 
% h_legend_r4 = legend(Legend_r4, 'interpreter','latex','location','best'); 
% 
% xlabel('$r$ [-]','Fontsize',14)
% ylabel('$h$ [-]','Fontsize',14)
% set(gca,'FontSize',14)
% 
% set(gcf, 'unit', 'inches');
% figure_size_4 =  get(gcf, 'position');
% 
% set(h_legend_r4, 'location', 'northeastoutside');
% set(h_legend_r4, 'unit', 'inches');
% legend_size = get(h_legend_r4, 'position');
% figure_size_4(3) = figure_size_4(3) + legend_size(3);
% set(gcf, 'position', figure_size_4);
% 
% set(hfig2b,'Units','Inches');
% pos = get(hfig2b,'Position');
% set(hfig2b,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
% print(hfig2b,'lastRuptureEvent_bothScaled','-dpdf','-r300');

%%

% hfig2c = figure;
% hfig2c.Renderer = 'Painters';
% 
% for j = 2:max(size(Out2))  % max(size(Out2))
%     newFolder = Out2{j};
%     cd(newFolder);
%     load('lastFrame.mat');
% %     all_figs = findobj(0, 'type', 'figure');
% %     delete(setdiff(all_figs, figs2keep));
%     plot(x*l_scale/R_film_final, h_store(:,end)*h0_init*10^9/(R_film_final*10^6), 'color', ColOrd(j,:), 'Linewidth', 1.5);
%     xlim([0 1])
%     Legend_r4{j} = strcat('$R_f = $', num2str(R_f*10^6),' $\mu$m');
%     hold on
%     cd ..
% end
% 
% h_legend_r4 = legend(Legend_r4{2:max(size(Out2))}, 'interpreter','latex','location','best'); 
% 
% xlabel('$r$ [-]','Fontsize',14)
% ylabel('$h$ [-]','Fontsize',14)
% set(gca,'FontSize',14)
% 
% set(gcf, 'unit', 'inches');
% figure_size_4 =  get(gcf, 'position');
% 
% set(h_legend_r4, 'location', 'northeastoutside');
% set(h_legend_r4, 'unit', 'inches');
% legend_size = get(h_legend_r4, 'position');
% figure_size_4(3) = figure_size_4(3) + legend_size(3);
% set(gcf, 'position', figure_size_4);
% 
% set(hfig2c,'Units','Inches');
% pos = get(hfig2c,'Position');
% set(hfig2c,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
% print(hfig2c,'lastRuptureEvent_bothScaled_R_film_final','-dpdf','-r300');

%%

% hfig2d = figure;
% hfig2d.Renderer = 'Painters';
% 
% for j = 2:max(size(Out2))  % max(size(Out2))
%     newFolder = Out2{j};
%     cd(newFolder);
%     load('lastFrame.mat');
% %     all_figs = findobj(0, 'type', 'figure');
% %     delete(setdiff(all_figs, figs2keep));
%     x = [0; x];
%     h_ruptureProf = [h_store(1,end); h_store(:,end)];
%     plot(x*l_scale/R_film_final, h_ruptureProf*h0_init*10^9/(2.734*R_film_final*10^6 - 31.08), 'color', ColOrd(j,:), 'Linewidth', 1.5);
%     xlim([0 1.2])
%     ylim([0 1.2])
%     Legend_r4{j} = strcat('$R_f = $', num2str(R_f*10^6),' $\mu$m');
%     hold on
%     cd ..
% end
% 
% h_legend_r4 = legend(Legend_r4{2:max(size(Out2))}, 'interpreter','latex','location','best'); 
% 
% xlabel('$r$ [-]','Fontsize',14)
% ylabel('$h$ [-]','Fontsize',14)
% set(gca,'FontSize',14)
% 
% set(gcf, 'unit', 'inches');
% figure_size_4 =  get(gcf, 'position');
% 
% set(h_legend_r4, 'location', 'northeastoutside');
% set(h_legend_r4, 'unit', 'inches');
% legend_size = get(h_legend_r4, 'position');
% figure_size_4(3) = figure_size_4(3) + legend_size(3);
% set(gcf, 'position', figure_size_4);
% 
% set(hfig2d,'Units','Inches');
% pos = get(hfig2d,'Position');
% set(hfig2d,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
% print(hfig2d,'ruptureProfile_finalR_film_scaling03','-dpdf','-r300');

%%

% hfig2e = figure;
% hfig2e.Renderer = 'Painters';
% 
% for j = 2:max(size(Out2))  % max(size(Out2))
%     newFolder = Out2{j};
%     cd(newFolder);
%     load('lastFrame.mat');
% %     all_figs = findobj(0, 'type', 'figure');
% %     delete(setdiff(all_figs, figs2keep));
%     x = [0; x];
%     h_ruptureProf = [h_store(1,end); h_store(:,end)];
%     plot(x*l_scale/R_film_final, h_ruptureProf/h_centre_j(end), 'color', ColOrd(j,:), 'Linewidth', 1.5);
%     xlim([0 1.0])
% %     ylim([0 1.2])
%     Legend_r4{j} = strcat('$R_f = $', num2str(R_f*10^6),' $\mu$m');
%     hold on
%     cd ..
% end
% 
% h_legend_r4 = legend(Legend_r4{2:max(size(Out2))}, 'interpreter','latex','location','best'); 
% 
% xlabel('$r$ [-]','Fontsize',14)
% ylabel('$h$ [-]','Fontsize',14)
% set(gca,'FontSize',14)
% 
% set(gcf, 'unit', 'inches');
% figure_size_4 =  get(gcf, 'position');
% 
% set(h_legend_r4, 'location', 'northeastoutside');
% set(h_legend_r4, 'unit', 'inches');
% legend_size = get(h_legend_r4, 'position');
% figure_size_4(3) = figure_size_4(3) + legend_size(3);
% set(gcf, 'position', figure_size_4);
% 
% set(hfig2e,'Units','Inches');
% pos = get(hfig2e,'Position');
% set(hfig2e,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
% print(hfig2e,'ruptureProfile_scales_h_centre_R_film_final','-dpdf','-r300');
% 
% 
% %%
% 
% hfig6a = figure;
% hfig6a.Renderer = 'Painters';
% 
% for j = 1:max(size(Out2))  % max(size(Out2))
%     newFolder = Out2{j};
%     cd(newFolder);
%     load('lastFrame.mat');
% %     all_figs = findobj(0, 'type', 'figure');
% %     delete(setdiff(all_figs, figs2keep));
%     x = [0; x];
%     min_Pres(j) = min(Pressure_rupture_m1);
%     h_ruptureProf = [h_store(1,end); h_store(:,end)];
%     plot(x(2:end-1)*l_scale/R_film_final, CapillaryPres_rupture_m1/abs(CapillaryPres_rupture_m1), 'color', ColOrd(j,:), 'Linewidth', 1.5);
%     xlim([0 2])
% %     ylim([0 1.2])
%     Legend_r4{j} = strcat('$R_f = $', num2str(R_f*10^6),' $\mu$m');
%     hold on
%     cd ..
% end
% 
% h_legend_r4 = legend(Legend_r4{1:max(size(Out2))}, 'interpreter','latex','location','best'); 
% 
% xlabel('$r$ [-]','Fontsize',14)
% ylabel('$C.P.$ [-]','Fontsize',14)
% set(gca,'FontSize',14)
% 
% set(gcf, 'unit', 'inches');
% figure_size_4 =  get(gcf, 'position');
% 
% set(h_legend_r4, 'location', 'northeastoutside');
% set(h_legend_r4, 'unit', 'inches');
% legend_size = get(h_legend_r4, 'position');
% figure_size_4(3) = figure_size_4(3) + legend_size(3);
% set(gcf, 'position', figure_size_4);
% 
% set(hfig6a,'Units','Inches');
% pos = get(hfig6a,'Position');
% set(hfig6a,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
% print(hfig6a,'CapPressureProfile_m1_scaled','-dpdf','-r300');
% %%
% 
% hfig7 = figure;
% hfig7.Renderer = 'Painters';
% 
% loglog(R_film_imposed, min_Pres, 'o')
% 
% xlabel('$R_{film}(t=0)$ [-]','Fontsize',14)
% ylabel('$minimum Pressure$ [-]','Fontsize',14)
% set(gca,'FontSize',14)
% 
% set(hfig7,'Units','Inches');
% pos = get(hfig7,'Position');
% set(hfig7,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
% print(hfig7,'minPress_Vs_R_film','-dpdf','-r300');

%%

hfig7 = figure;
hfig7.Renderer = 'Painters';

for j = 1:max(size(Out2))  % max(size(Out2))
    newFolder = Out2{j};
    cd(newFolder);
    load('lastFrame.mat');
    x = [0; x];
    min_Pres(j) = min(Pressure_rupture_m1);
    h_ruptureProf = [h_store(1,end); h_store(:,end)];
    plot(x(2:end-1)*l_scale/R_film_final, Pres_h_max, 'color', ColOrd(j,:), 'Linewidth', 1.5);
    xlim([0 2])
    Legend_r4{j} = strcat('$R_f = $', num2str(R_f*10^6),' $\mu$m');
    hold on
    cd ..
end

toc
