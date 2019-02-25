clear all
close all
clc

format long g
tic

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0, 'defaulttextInterpreter','latex');

filmConfiguration = 'flatFilms_PBC';
% filmConfiguration = 'semiInfiniteNonFlatFilms';
% filmConfiguration = 'finiteSizedNonFlatFilms';

switch filmConfiguration
    case 'flatFilms_PBC' 
        kappa = 0;
        transitionLength = 0;
        L_flat = 15;
        L_curv = 0;
     
        deltaX = 0.05;
        N = (L_flat+L_curv)./deltaX;
        ctimestep = 2.75;
        deltaT = deltaX^ctimestep;
        endTime = 150;
        seN = 20;
        N_Reals = 2;                        % number of realizations
        startRealization = 1;

        
        Tmp = 0.001;
        [h x] = initialProfile(kappa,L_flat,L_curv,transitionLength,deltaX, filmConfiguration);
    case 'semiInfiniteNonFlatFilms'
        kappa = 0.1;
        transitionLength = 1./sqrt(2*kappa);
        L_flat = 240;
        L_curv = 105;
        
        deltaX = 0.05;
        N = (L_flat+L_curv)./deltaX;
        ctimestep = 2.75;
        deltaT = deltaX^ctimestep;
        endTime = 150;
        seN = 20;
        N_Reals = 1;                        % number of realizations
        startRealization = 1;


        Tmp = 0.00;
        [h x] = initialProfile(kappa,L_flat,L_curv,transitionLength,deltaX, filmConfiguration);
    case 'finiteSizedNonFlatFilms'
        %% raw conditions from the paper
        R_film = [20 25 30 35 40 50 60 65 70 75 80 85 90 100 115 150 200 300 400 500 600 700 800 900 1000];   % radius of the film
        % R_film = [30];
        R_f = R_film.*10^-6;              % in m
        h0_init = 150e-9;                 % initial film height in m
        A_vw = 1.25e-21;                  % Hamaker constant
        gam = 0.034;                      % surface tension
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
        hJoyeStart = 1;                                 % determine when to start measuring thinning rates
        hJoyeEnd = 0.627*kappa^(-2/7);                  % where to end
        h_critical_start = 0.627*kappa^(-2/7)*1.2;      % this is when the film thinning velocity starts becoming nearly constant
        h_critical_end = 0.627*kappa^(-2/7)*0.8;        % this is when the film thinning velocity ends becoming nearly constant


        %% domain size and discretization parameters

        deltaX = 0.00125*ones(size(R_f));            % grid size (tested for grid independent results)
        L_flat = round(R_f/l_scale,4);              % length of the flat film
        N_flat = round(L_flat./deltaX);             % number of grid points in the same
        ctimestep = 2.75;                              % exponent used in deciding deltaT = deltaX^c --> although c = 2.75 suffices, but a higher temp resolution enables more time stamps

        for i = 1:length(N_flat)
            if (N_flat(i)) < 40
                N_flat(i) = 40;
                deltaX(i) = L_flat(i)/N_flat(i);
                deltaT(i) = deltaX(i)^ctimestep;                      % end time of a realization
            end        
        end
        deltaT = deltaX.^ctimestep;


        %% simulation parameters
        Tmp = 0.0;                          % dimensionless noise strength (= 0, for deterministic)
        upperLimitOnL_curv = sqrt(pi*h0_init^2*gam/(2*A_vw*kappa));
        lowerLimitOnL_curv = 1./sqrt(2*kappa);
        transitionLength = lowerLimitOnL_curv;
        L_curv = 0.1;                       % length of the curved portion of the film, for kappa > 1, one needs a smaller value of of L_curv
        endTime = 150;
        seN = 1;                            % save every these many time steps
        N_Reals = 1;                        % number of realizations
        animationSkip = 20;                 % save animation every these many time steps
        startRealization = 1;               % first realization

        for i = 1:length(L_flat)
            if L_flat(i) <= transitionLength
                error('Length of the flat film is smaller than the transition region')
                break;
            end
        end
        
        %% plot the first film length
        [h x] = initialProfile(kappa,L_flat(1),L_curv,transitionLength,deltaX(1),filmConfiguration);
        plot(x,h,'o')
        ylim([0.95 5])
        
        %% parent folder
        
        mk = strcat('h0_',num2str(h0_init*10^9),'nm','_Avw_',num2str(A_vw),'_ST_',num2str(gam),'_Rc_',num2str(Rc));
        mkdir(mk);
        cd(mk);
        
        %% start simulations for different radii of films

        for i = 1:length(L_flat)
            str1{i} = strcat('Rf_', num2str(R_film(i)),'_mu_m');
            mkdir(str1{i});
            path_dest{i} = strcat('./', num2str(str1{i}));
            destinatn = path_dest{i};
            copyfile('../*.m', destinatn)
            run_mainFiles{i} = strcat('./',num2str(str1{i}));
            cd(run_mainFiles{i})
            main_finiteSizedFilms(filmConfiguration, R_f(i), h0_init, A_vw, gam, Rc, visc, L_flat(i), N_flat(i), deltaX(i), deltaT(i), transitionLength, h_drain_start, h_drain_end, h_critical_start,...
                h_critical_end, t_cr_dimensional, res_limit, ctimestep, Tmp, L_curv, endTime, seN,...
                N_Reals, animationSkip, startRealization, hJoyeStart, hJoyeEnd)
            
            cd ..
            fileToBeSaved = strcat('workspace_','Rf_',num2str(min(R_film)),'_to_',num2str(max(R_film)),'.mat');
            save('fileToBeSaved')
        end
end
        
%% start simulations for semi-infinite films or flat films with periodic boundary conditions

for i = 1:length(L_flat)
    switch filmConfiguration
        case {'flatFilms_PBC','semiInfiniteNonFlatFilms'}
            str1{i} = strcat('theta_', num2str(Tmp));
            mkdir(str1{i});
            cd(str1{1});
            str2{i} = strcat('kappa_', num2str(kappa));
            mkdir(str2{i});
            path_dest{i} = strcat('./', num2str(str2{i}));
            destinatn = path_dest{i};
            copyfile('../*.m', destinatn)
            run_mainFiles{i} = strcat('./',num2str(str2{i}));
            cd(run_mainFiles{i})
            main_flatOrSemiInfiniteFilms(filmConfiguration, kappa, L_flat, deltaX, deltaT, transitionLength, ctimestep, Tmp, L_curv,...
                endTime, seN, N_Reals, startRealization)
            cd ..
    end
end

save('workspaceInputParameters.mat')

cd ..

toc

%%
%%%%%%%%%%%%%%%%%%%%%% end this file here  %%%%%%%%%%%%%%%%%%%%%%%%%%%




