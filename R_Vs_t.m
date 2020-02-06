clear all
close all
clc

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaulttextInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');

tic

relFolder = {'axisSymmetricFilm\disjPress_on\repulsion_off\Unadj_CP_h0_300nm_Avw_1.5e-20_ST_0.0445_Rc_0.0018_200_50nm_c1_0_c2_0_June13th', ...
             'axisSymmetricFilm\disjPress_on\repulsion_off\Unadj_CP_h0_500nm_Avw_1.5e-20_ST_0.0445_Rc_0.0018_200_50nm_c1_0_c2_0_June13th', ...
             'axisSymmetricFilm\disjPress_on\repulsion_off\Unadj_CP_h0_1000nm_Avw_1.5e-20_ST_0.0445_Rc_0.0018_200_50nm_c1_0_c2_0_June13th', ...
             'axisSymmetricFilm\disjPress_on\repulsion_off\Unadj_CP_h0_2000nm_Avw_1.5e-20_ST_0.0445_Rc_0.0018_200_50nm_c1_0_c2_0_June13th'}; 

h0 = [300 500 1000 2000];

for j = 1:length(h0)
%     pathFol{j} = relFolder{j};
    cd(relFolder{j})
    load('results_differentFilmSize_diffThinnRates.mat')
    R_film_tr{j} = R_film_final*10^6;
    R_film_t0{j} = R_film;
    ratio_R_film{j} = R_film_final*10^6./R_film;
    cd('../../../..')
end
    
% R_t0 = [30 50 75 100 150 200 300 500 1000];
% R_t0_2000 = [50 75 100 150 200 300 500 1000];
% R_t_by_R_t0 = [25 48 74 99 149 199 299 499 999; 23 46 72 99 149 199 299 499 999; 16 42 71 96 147 197 298 499 999]./R_t0;
% R_t_by_R_t0_2000 = [35 63 93 144 197 299 499 999]./R_t0_2000;

h1a = figure;
h1a.Renderer = 'Painters';

% semilogx(R_t0, R_t_by_R_t0, 'o')
% hold on
% semilogx(R_t0_2000, R_t_by_R_t0_2000, 'o')
for i = 1:length(R_film_tr)
    semilogx(R_film_t0{i}, ratio_R_film{i},'o');
    hold on
end
ylim([0.05 1.1]);
xlim([10 4100]);
xlabel('$R_{film}$ ($t = t_0$)($\mu$m)','Fontsize',14);
ylabel('$R_{film}$($t = t_r$)/$R_{film}$($t = t_0$) [-]','Fontsize',14);
set(gca,'FontSize',14);
legend('$h_o = 300$ nm', '$h_o = 500$ nm', '$h_o = 1000$ nm', '$h_o = 2000$ nm', 'Location', 'best');

set(h1a,'Units','Inches');
pos = get(h1a,'Position');
set(h1a,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(h1a,'ratioOfFinalRadiustoInitialRadius_04','-dpdf','-r300');

toc