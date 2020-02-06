function [t_rupt drainageTime drainageTime_right drainageTime_right_rupt ...
    avg_cr_thinningRate_fit h_cr_final h_cr_final_FullFilmavg] = main_finiteSizedFilms(filmConfiguration, disjPress_switch, R_f, h0_init, A_vw, gam, Rc, visc, L_flat, N_flat, deltaX, deltaT, transitionLength, h_drain_start, h_drain_end, h_critical_start,...
                h_critical_end, t_cr_dimensional, res_limit, ctimestep, Tmp, L_curv, endTime, seN, N_Reals, animationSkip, ...
                startRealization, hJoyeStart, hJoyeEnd, cutOff_thickness, correctionLP_switch, repulsion_coeff, repulsion_expon, vdW_repulsion,...
                eq_thickness_EDL_vdW)



%% derived quantities (have also been derived in the input file)
kappa = pi*h0_init^3*gam/A_vw/Rc;
t_scale = 12*pi^2*visc*gam*h0_init^5/A_vw^2;
l_scale = h0_init^2*sqrt(2*pi*gam/A_vw); 
t_cr = t_cr_dimensional./t_scale;

%% simulation set-up

[h x] = initialProfile(kappa,L_flat,L_curv, R_f, Rc,transitionLength,deltaX, filmConfiguration, correctionLP_switch);  
length(h)
L = 2*(L_curv + L_flat);   % total length of the film (curved+flat)
N = length(x) - 1;   %  -1 (to keep the notion consistent with N being the number of intervals and not the number of grid points);  
gx = gx_generator(N,L,x);  % generates a matrix that is going to be used when we finally implement noise


hfig = figure;
subplot(2,1,1)
plot(x,h)
% ylim([0 5])
subplot(2,1,2)
plot(x,h,'o')
ylim([0.95 1.5])

set(hfig,'Units','Inches');
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hfig, 'initialFilmProfile' ,'-dpdf','-r300')


%% to check symmetry in the initial condition
% h_left = h(1:(length(h)-1)/2);
% h_right = flipud(h((length(h)+3)/2:end));
% h_diff = h_left - h_right;

%%
realization = (startRealization-1);                % counter for the number of realizations
t_rupt = zeros(N_Reals,1);      % preallocate rupture times vector

for m = 1:N_Reals      
    h_adjusted= N + 1; % extra 1 for the conversion from the number of intervals to number of grid points (2 ghost points are already included in the initial profile)
    A=spalloc( h_adjusted , h_adjusted ,10 + h_adjusted *5-10) ;% preallocate sparse matrix
    %% generate a 'p' vector that is going to be used to fill up band in the sparse matrix
    p1=linspace (3 , h_adjusted-2,h_adjusted-4) *(1 + h_adjusted ) ;
    p = [] ;
    k = linspace(1,h_adjusted,h_adjusted);
    for i =-3:1
        p=[p p1+ i*h_adjusted ] ;
    end
    %% boundary conditions : dh/dx = 0 and d3h/dx3 = 0 for right far field || and h = 1 + kappa*x^2 and d2h/dx2 = 2 kappa for the far left (curved) film
    p_ones=[1 2*h_adjusted+1 h_adjusted+2 h_adjusted^2];
    A(p_ones )=ones(1,4);
    %p_minusOnes=[h_adjusted*(h_adjusted-4) (h_adjusted-1)*(h_adjusted)+(h_adjusted-1)];
    %A(p_minusOnes)=ones (1,2) *-1;
    A(h_adjusted+1) = -2;
    %A(h_adjusted*(h_adjusted-3)) = 3;
    %A(h_adjusted*(h_adjusted-1)) = -3;
    A(h_adjusted*(h_adjusted-2)+(h_adjusted-1))=1;
    A(h_adjusted*(h_adjusted-2))=1;
    A(h_adjusted*(h_adjusted-1))=-2;
    %% call the solver

    [t_rupt(m) x_rupt(m) minH(:,m)] = filmSolver(filmConfiguration, disjPress_switch,correctionLP_switch, repulsion_coeff, repulsion_expon, vdW_repulsion,...
                                    L_flat, R_f, Rc, transitionLength,L_curv,N,deltaX,deltaT,kappa,Tmp,gx,h_adjusted,A,p,endTime,seN, res_limit, cutOff_thickness, eq_thickness_EDL_vdW);

    reali_series(m) = m;
    realization = realization + 1;
    mk = strcat('realization',num2str(realization));  % name your realization folder
    mk2 = mkdir(mk);                                  % make its directory
    movefile('hData.mat',mk)                              % move all the data files to that directory
end

    
%% store the data (but more importantly the rupture times) into a mat file, so that there is no further post processing required if we are just looking for T_r
filename = ['data_', 'kappa_',num2str(kappa),'_Lf_',num2str(L),'_N_',num2str(N), '_Tmp_', num2str(Tmp),'.mat'];


[t_rupt drainageTime drainageTime_right drainageTime_right_rupt  ...
    avg_cr_thinningRate_fit h_cr_final h_cr_final_FullFilmavg] = postProc_det(filmConfiguration, disjPress_switch, correctionLP_switch, R_f, L_flat, L_curv, transitionLength, ...
    deltaX, deltaT, kappa, t_scale, l_scale, seN, animationSkip, h_drain_start, h_drain_end,...
    h_critical_start, h_critical_end, t_cr, res_limit, hJoyeStart, hJoyeEnd, h0_init, Rc, gam);

save(filename)


end
