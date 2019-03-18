clear all
close all
clc

format long g
tic

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0, 'defaulttextInterpreter','latex');

%% Select the film configuration that you want to simulate

%%%% flat films with periodic boundary conditions %%%%%%
% filmConfiguration = 'flatFilms_PBC';    

%%%%  semi-infinite films with far field boundary conditions %%%%%%
filmConfiguration = 'semiInfiniteNonFlatFilms';

%%%%% finite sized films with curvature on both sides  %%%%%%
% filmConfiguration = 'finiteSizedNonFlatFilms';

%% Do we switch off the disjoining pressure in the solver?

disjPress_switch = 'on';


%% Based on the choice a part of this code gets executed

switch filmConfiguration
    case {'flatFilms_PBC', 'semiInfiniteNonFlatFilms'} 
        kappa = 0.001;
        Tmp = 0.0;
        L_flat = 240;  
        if kappa == 0
            transitionLength = 0;
            L_curv = 0;
        else
            transitionLength = 1./sqrt(2*kappa);
            L_curv = 300;
        end

        if Tmp == 0
            N_Reals = 3;
        else
            N_Reals = 2;   %% please adjust the number of realizations based on how many you want to sample
        end

        if isequal(disjPress_switch, 'on') 
            cutOff_thickness = 1e-05;            % keep it lower to be able to probe even smaller thicknesses if it reaches 
        elseif isequal(disjPress_switch, 'off') 
            cutOff_thickness = 0.05;             % since the thinning rate in the absence of disj pres decreases asymptotically, a higher cut-off would save computational time
        end
     
        deltaX = 0.05;
        N = (L_flat+L_curv)./deltaX;
        ctimestep = 2.75;
        deltaT = deltaX^ctimestep;
        endTime = 150;
        seN = 20;
        startRealization = 1;

        [h x] = initialProfile(kappa,L_flat,L_curv,transitionLength,deltaX, filmConfiguration);

        %% here we kick-start the simulations

        str1 = strcat('theta_', num2str(Tmp));
        mkdir(str1);
        cd(str1);
        str2 = strcat('kappa_', num2str(kappa));
        mkdir(str2);
        path_dest = strcat('./', num2str(str2));
        destinatn = path_dest;
        copyfile('../*.m', destinatn)
        run_mainFiles = strcat('./',num2str(str2));
        cd(run_mainFiles)
        main_flatOrSemiInfiniteFilms(filmConfiguration,disjPress_switch, kappa, L_flat, deltaX, deltaT, transitionLength, ctimestep, Tmp, L_curv,...
            endTime, seN, N_Reals, startRealization, cutOff_thickness)
        cd ..
    case 'finiteSizedNonFlatFilms'
        %% raw conditions from the paper
        R_film = [25 30 35 40 50 60 65 70 75 80 85 90 100 115 150 200 300 400 500 600 700 800 900 1000];   % radius of the film
%         R_film = [70 90 150 800];
        R_f = R_film.*10^-6;              % in m
        h0_init = 1000e-9;                % initial film height in m
        A_vw = 1.5e-20;                   % Hamaker constant
        gam = 0.0445;                     % surface tension
        Rc = 1.8e-3;                      % radius of capillary
        visc = 0.00089;                   % viscosity
        t_cr_dimensional = [2 3 4 5];     % time resolution (in sec) 
        res_limit_dimensional = 22e-6;    % spatial resolution (in m)

        %% derived quantities from above 

        l_scale = h0_init^2*sqrt(2*pi*gam/A_vw);        % length scale of the system - obtained from the O(1) scaling
        kappa = pi*h0_init^3*gam/A_vw/Rc;               % dimensionless curvature - the free parameter of the system
        t_scale = 12*pi^2*visc*gam*h0_init^5/A_vw^2;    % time scale of the system

        %% drainage time start and end, critical thinning rates start and end, Joye's thinning rate start and end

        res_limit = res_limit_dimensional/l_scale;      % spatial resolution (scaled)
        h_drain_start = 100e-9/h0_init;                 % 100nm as mentioned in Wasan & Malhotra (1987)
        h_drain_end = 25e-9/h0_init;                    % 25 nm as mentioned in Wasan & Malhotra (1987)
        hJoyeStart = 0.7;                                 % determine when to start measuring thinning rates
        hJoyeEnd = 0.627*kappa^(-2/7);                  % where to end
        h_critical_start = 0.627*kappa^(-2/7)*1.2;      % this is when the film thinning velocity starts becoming nearly constant
        h_critical_end = 0.627*kappa^(-2/7)*0.8;        % this is when the film thinning velocity ends becoming nearly constant


        %% domain size and discretization parameters

        deltaX = 0.0025*ones(size(R_f));            % grid size (tested for grid independent results)
        L_flat = round(R_f/l_scale,4);              % length of the flat film
        N_flat = round(L_flat./deltaX);             % number of grid points in the same
%         ctimestep = 3;                              % exponent used in deciding deltaT = deltaX^c --> although c = 2.75 suffices, but a higher temp resolution enables more time stamps
        for i = 1:length(R_film)
            if R_film(i) < 50
                seN(i) = 20;                            % save every these many time steps
                ctimestep(i) = 3;
            else
                seN(i) = 20;                            % save every these many time steps
                ctimestep(i) = 3.0;
            end
        end

        for i = 1:length(N_flat)
            if (N_flat(i)) < 40
                N_flat(i) = 40;
                deltaX(i) = L_flat(i)/N_flat(i);
                deltaT(i) = deltaX(i).^ctimestep(i);     % end time of a realization
            end        
        end
        deltaT = deltaX.^ctimestep;


        %% simulation parameters
        
        Tmp = 0.0;                          % dimensionless noise strength (= 0, for deterministic)
        upperLimitOnL_curv = sqrt(pi*h0_init^2*gam/(2*A_vw*kappa));
        lowerLimitOnL_curv = 1./sqrt(2*kappa);
        transitionLength = 0.5*lowerLimitOnL_curv;
        L_curv = 15*lowerLimitOnL_curv;                       % length of the curved portion of the film, for kappa > 1, one needs a smaller value of of L_curv
        endTime = 0.0001;
        
        N_Reals = 1;                        % number of realizations
        animationSkip = 50;                 % save animation every these many time steps
        startRealization = 1;               % first realization

        for i = 1:length(L_flat)
            if L_flat(i) <= transitionLength
                error('Length of the flat film is smaller than the transition region')
                break;
            end
        end
        
        if isequal(disjPress_switch, 'on') 
            cutOff_thickness = 1e-05;            % keep it lower to be able to probe even smaller thicknesses if it reaches 
        elseif isequal(disjPress_switch, 'off') 
            cutOff_thickness = 0.05;             % since the thinning rate in the absence of disj pres decreases asymptotically, a higher cut-off would save computational time
        end
        
        %% plot the first film length
        [h x] = initialProfile(kappa,L_flat(1),L_curv,transitionLength,deltaX(1),filmConfiguration);
        plot(x,h,'o')
        ylim([0.95 5])
        
        %% parent folder
        
        mk = strcat('h0_',num2str(h0_init*10^9),'nm','_Avw_',num2str(A_vw),'_ST_',num2str(gam),'_Rc_',num2str(Rc), '_disjPr_',disjPress_switch);
        mkdir(mk);
        cd(mk);
        copyfile('../*.m', '.')
        copyfile('../dataRadoev1984.xlsx', '.')
        
        %% start simulations for different radii of films

        for i = 1:length(L_flat)
            str1{i} = strcat('Rf_', num2str(R_film(i)),'_mu_m');
            mkdir(str1{i});
            path_dest{i} = strcat('./', num2str(str1{i}));
            destinatn = path_dest{i};
            copyfile('*.m', destinatn)
            run_mainFiles{i} = strcat('./',num2str(str1{i}));
            cd(run_mainFiles{i})
            [t_rupt(i) drainageTime(i) drainageTime_left(i) drainageTime_right(i) drainageTime_right_rupt(i) drainageTime_left_rupt(i) ...
            avg_cr_thinningRate_fit(i) h_cr_final(:,i) h_cr_final_FullFilmavg(:,i)] = main_finiteSizedFilms(filmConfiguration, disjPress_switch, ...
                        R_f(i), h0_init, A_vw, gam, Rc, visc, L_flat(i), N_flat(i), deltaX(i), deltaT(i), transitionLength, h_drain_start,...
                        h_drain_end, h_critical_start, h_critical_end, t_cr_dimensional, res_limit, ctimestep(i), Tmp, L_curv, endTime, seN(i),...
                        N_Reals, animationSkip, startRealization, hJoyeStart, hJoyeEnd, cutOff_thickness);

            
            cd ..
            fileToBeSaved = strcat('workspace_','Rf_',num2str(min(R_film)),'_to_',num2str(max(R_film)),'.mat');
            save(fileToBeSaved)
        end
        save('results_differentFilmSize.mat')
        
        %% post process finite sized films
        
        postProcessFiniteSizedFilms();
end

save('workspaceInputParameters.mat')

cd ..

toc

%%
%%%%%%%%%%%%%%%%%%%%%% end this file here  %%%%%%%%%%%%%%%%%%%%%%%%%%%




