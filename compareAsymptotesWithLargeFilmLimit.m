clear all
close all
clc

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaulttextInterpreter','latex');

h0 = [3e-7 5e-7 1e-6 2e-6];
gam = 0.0445;
A_vw = 1.5e-20;
Rc = 1.8e-3;
visc = 0.00089;
kappa = pi*h0.^3*gam./(A_vw*Rc);

%% relevant time scales 

t_scale_smallFilm = 3/2*visc/gam*(Rc^2./h0);
t_scale_largeFilm_semiInfFilms = 12*pi^2*gam*visc*h0.^5/A_vw^2;

%% dimensionless rupture times from axisymmetric simulations
t_r = [22.87 35.74 71.28 150.48];
% t_r = [20.46 70.84 184.065 730.13 3058.8];

% making it dimensional
t_r_dimensionless = t_r./t_scale_smallFilm;

t_r_model = 6.32.*kappa.^(4/7)

hfig10 = figure;
hfig10.Renderer = 'Painters';

loglog(kappa, t_r_model, '-')
hold on
loglog(kappa, t_r_dimensionless, 'o')
ylim([50 5000])
xlabel('$\kappa$ (-)','Fontsize',18)
ylabel('$\tilde{t}_r$ (-)','Fontsize',18)
set(gca,'FontSize',18)

set(hfig10,'Units','Inches');
pos = get(hfig10,'Position');
set(hfig10,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hfig10,'fig10','-dpdf','-r300')

%%



% converting it to equivalent dimensionless rupture times for semi-infinite
% films
t_r_LargeFilmLimit = t_r./t_scale_largeFilm_semiInfFilms;

%% dimensional rupture times if \tilde{t}_r = kappa^(-10/7) for axissym case

t_r_reAdjusted = 432*pi^(4/7)*visc*gam^(-3/7).*h0.^(5/7)./(A_vw^(4/7)*Rc^(-10/7));

%% dimensional rupture times if \tilde{t}_r = kappa^(-10/7) for 2D case

t_r_adjusted = 1.05*12*pi^(4/7)*visc*gam^(-3/7).*h0.^(5/7)./(A_vw^(4/7)*Rc^(-10/7));

%% dimensionless rupture times model from PRF paper
t_r_largeFilm = 3.92./(1 + 3.74*kappa.^(10/7));
% t_r_largeFilm_sim = [0.000596 0.000069499947 0.000003375929 7.1859e-08];

% ratio_dimensional = t_r_dimensional./t_r_adjusted
% ratio_dimensionless = t_r_LargeFilmLimit./t_r_largeFilm


hfig0 = figure;
hfig0.Renderer = 'Painters';

loglog([1 1 1 1], t_r_adjusted, '+')
ylim([9.6795 1.1184e3])
xlabel('$h_o$ [-]','Fontsize',18)
ylabel('$t_r$ [-]','Fontsize',18)
set(gca,'FontSize',18)

set(hfig0,'Units','Inches');
pos = get(hfig0,'Position');
set(hfig0,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hfig0,'largeFilmLimit','-dpdf','-r300')

hfig = figure;
hfig.Renderer = 'Painters';

loglog(kappa, t_r_dimensionless, 'o')
% hold on
% loglog(kappa, t_r_largeFilm_sim, '+')
hold on
loglog(kappa, t_r_reAdjusted)
xlabel('$h_o$ [-]','Fontsize',18)
ylabel('$t_r$ [-]','Fontsize',18)
set(gca,'FontSize',18)

set(hfig,'Units','Inches');
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hfig,'compareLargeFilmLimits_dimensional','-dpdf','-r300')

hfig2 = figure;
hfig2.Renderer = 'Painters';

loglog(kappa, t_r_LargeFilmLimit, 'o')
% hold on
% loglog(kappa, t_r_largeFilm_sim, '+')
% hold on
% loglog(kappa, t_r_largeFilm)
hold on
loglog(kappa, t_r_largeFilm)
xlabel('$\kappa$ [-]','Fontsize',18)
ylabel('$\tilde{t}_r$ [-]','Fontsize',18)
set(gca,'FontSize',18)

set(hfig2,'Units','Inches');
pos = get(hfig2,'Position');
set(hfig2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hfig2,'compareLargeFilmLimits_dimensionless','-dpdf','-r300')



