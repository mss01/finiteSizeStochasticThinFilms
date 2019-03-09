function [v_re_det t_re v_MTR t_MTR v_MTR_1997Paper v_MTR_Tsekov] = Reynolds_and_MTR(h_dimensionless, kappa, L_flat, R_f, h0_init,...
                                                                    t_scale, h_drain_start, h_drain_end, visc, gam, Rc, A_vw)

v_re_det = ((1 + 6.*h_dimensionless.^3*kappa)./(L_flat).^2)'*h0_init/t_scale*10^10;
preFactorFunc = @(x) 1./(1+6.*kappa.*x.^3);
preFactorRe = integral(preFactorFunc, h_drain_end ,h_drain_start);
t_re = preFactorRe.*L_flat.^2.*t_scale;             %% preFactorRe is coming from integrating 1/(1 + 6h^3*kappa) from 100nm to 25nm

% MTR theory
% v_MTR = (((1 + 6*kappa)^8.*h_dimensionless.^12)./(3.*L_film_scaled.^4)).^0.2*h0_init/t_scale*10^10;
for i = 1:length(h_dimensionless)
    v_MTR(:,i) = 1/6.*(((1 + 6*kappa.*h_dimensionless(i).^3).^8)./(108.*h_dimensionless(i).^12*L_flat.^4)).^0.2*h0_init/t_scale*10^10;
end

preFactorFuncMTR = @(y) 6.*(108.*y.^12).^(1/5)./(1+6.*kappa.*y.^3).^(8/5);
% preFactorFuncMTR = @(y) 1./(1+6.*kappa.*y.^3);

preFactorMTR = integral(preFactorFuncMTR, h_drain_end ,h_drain_start);
t_MTR = preFactorMTR.*(L_flat.^4.).^(1/5).*t_scale;

h_MTR = h0_init.*h_dimensionless;
for i = 1:length(h_dimensionless)
        v_MTR_1997Paper(:,i) = 8/(3*visc)*((4.*h_MTR(i).^12.*(gam/Rc + A_vw/(6*pi*h_MTR(i).^3)).^8)./(3.83^12*gam^3.*(R_f).^4)).^0.2.*10^10;%% from Manev 1997 paper
        v_MTR_Tsekov(:,i) = 1/(6*visc).*((h_MTR(i).^12*(gam/Rc + A_vw/(6*pi*h_MTR(i).^3)).^8)./(4*gam^3.*(R_f).^4)).^0.2.*10^10;
end

end