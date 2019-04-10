function [t_rupt drainageTime drainageTime_right drainageTime_right_rupt  avg_cr_thinningRate_fit h_cr_final ...
                            h_cr_final_FullFilmavg] = postProc_det(filmConfiguration, disjPress_switch, R_f, L_flat, L_curv, transitionLength, ...
                                                                    deltaX, deltaT, kappa, t_scale, l_scale, seN, animationSkip, h_drain_start, h_drain_end, ...
                                                                    h_critical_start, h_critical_end, t_cr, res_limit, hJoyeStart, hJoyeEnd, h0_init, Rc);

tic

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0, 'defaulttextInterpreter','latex');

folders = dir('realization*');
load('realization1/hData.mat')
s = 1;
start = 0;
last = 0;

%% simulation set up

switch filmConfiguration
    case 'finiteSizedNonFlatFilms'
        L = 2*(L_curv + L_flat);                    % total length of the film (curved+flat)
    case 'axisSymmetricFilm'
        L = L_curv + L_flat;                        % total length of the film (curved+flat)
end
[h x] = initialProfile(kappa,L_flat,L_curv,transitionLength,deltaX, filmConfiguration);       % get the initial conditions
tt = seN*deltaT;                                            % the rate at which files were saved
cr_thickness = 0.627*kappa^(-2/7);                          % theoretical prediction of critical thickness  

%% identification of the extent of the flat films, usually where the dimple forms

[locDimple_right x_centre] = keyLocationsInFilm(filmConfiguration, x, L_flat, deltaX);
% [locDimple_left locDimple_right x_centre] = keyLocationsInFilm(x, L_flat, deltaX);

%% from here on starts the averaging over multiple realizations
[h_right_avg_j h_min h_centre_j v_thin_rim v_thin_centre avg_cr_thinningRate_fit h_cr_final ...
    h_cr_final_FullFilmavg drainageTime drainageTime_right drainageTime_right_rupt t_rupt h_max_dimp_r beginDrainageTime_right endDrainageTime_right x_dimple_loc_right] = ...
                                    loopingOverRealizations(x, L_flat, transitionLength, tt, x_centre, locDimple_right, res_limit, cr_thickness, ...
                                                                    h_drain_start, h_drain_end, t_cr, deltaX, deltaT, h_critical_start, h_critical_end);

% [h_right_avg_j h_left_avg_j h_avg h_min h_centre_j v_thin_rim v_thin_centre avg_cr_thinningRate_fit h_cr_final ...
%     h_cr_final_FullFilmavg drainageTime drainageTime_left drainageTime_right drainageTime_right_rupt drainageTime_left_rupt t_rupt h_max_dimp_l h_max_dimp_r] = ...
%                                     loopingOverRealizations(x, L_flat, transitionLength, tt, x_centre, locDimple_left, locDimple_right, res_limit, cr_thickness, ...
%                                                                     h_drain_start, h_drain_end, t_cr, deltaX, deltaT, h_critical_start, h_critical_end);
    
[h_min_rim h_centre_Joye t_rim v_re_Joye dhdt_rim dhdt_centre c_r_Joye ratio_v_vre ratio_vc_vre] = ...
    joyeAnalysis(filmConfiguration, disjPress_switch, hJoyeStart, hJoyeEnd, h_min, h_max_dimp_r, h_centre_j, deltaT, seN, t_store, kappa, L_flat, R_f, h0_init, Rc);
                        
save('workspace_deterministic_t_cr.mat')
% makeAnimation_det(filmConfiguration, animationSkip,kappa, L_flat, L_curv, transitionLength,deltaX, h_store, t_store, h0_init, t_scale, l_scale, beginDrainageTime_right, endDrainageTime_right, x_dimple_loc_right, res_limit);

toc

close all

end
