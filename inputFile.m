clear all
close all
clc

format long g
tic

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0, 'defaulttextInterpreter','latex');

%% Select the film configuration that you want to simulate

%%%%% flat films with periodic boundary conditions %%%%%%
% filmConfiguration = 'flatFilms_PBC';    

%%%%%  semi-infinite films with far field boundary conditions %%%%%%
% filmConfiguration = 'semiInfiniteNonFlatFilms';

%%%%% finite sized films with curvature on both sides  %%%%%%
% filmConfiguration = 'finiteSizedNonFlatFilms2D';

%%%%% axis-symmetric case  %%%%%%
filmConfiguration = 'axisSymmetricFilm';

%% Do we switch off the disjoining pressure in the solver?

disjPress_switch = 'on';
repulsion_switch = 'off';
adhocRepulsion_switch = 'off';

%% correction for Laplace Pressure (gets important when R_film becomes of the order of Rc)

correctionLP_switch = 'off';

%% switch for experimental or parametric study within axisSymmetric films

% ExpOrPar = 'parametric';
ExpOrPar = 'experimental';


%% Based on the choice a part of this code gets executed

switch filmConfiguration
    case {'flatFilms_PBC', 'semiInfiniteNonFlatFilms'} 
        mkdir(filmConfiguration)
        cd(filmConfiguration)
        kappa = 100;
        Tmp = 0.0;
        L_flat = 1; 
        R_f = 0;
        Rc = 0;
        repulsion_coeff = 0;
        repulsion_expon = 0;
        res_limit = 0;
        eq_thickness_EDL_vdW = 0;
        vdW_repulsion = 0;
        if kappa == 0
            transitionLength = 0;
            L_curv = 0;
        else
            transitionLength = 1./sqrt(2*kappa);
%             L_curv = 15*transitionLength;
            L_curv = 1;
        end

        if Tmp == 0
            N_Reals = 1;
        else
            N_Reals = 4;   %% please adjust the number of realizations based on how many you want to sample
        end

        if isequal(disjPress_switch, 'on') 
            cutOff_thickness = 1e-05;            % keep it lower to be able to probe even smaller thicknesses if it reaches 
            mkdir('disjPress_on')
            cd('disjPress_on')
        elseif isequal(disjPress_switch, 'off') 
            cutOff_thickness = 0.05;             % since the thinning rate in the absence of disj pres decreases asymptotically, a higher cut-off would save computational time
            mkdir('disjPress_off')
            cd('disjPress_off')
        end
     
        deltaX = 0.005;
        N = (L_flat+L_curv)./deltaX;
        ctimestep = 3.25;
        deltaT = deltaX^ctimestep;
        endTime = 0.01; %150;
        seN = 5;
        startRealization = 1;


        %% here we kick-start the simulations

        str1 = strcat('theta_', num2str(Tmp));
        mkdir(str1);
        cd(str1);
        str2 = strcat('kappa_', num2str(kappa));
        mkdir(str2);
        path_dest = strcat('./', num2str(str2));
        destinatn = path_dest;
        copyfile('../../../*.m', destinatn)
        run_mainFiles = strcat('./',num2str(str2));
        cd(run_mainFiles)
        main_flatOrSemiInfiniteFilms(filmConfiguration, disjPress_switch, kappa, L_flat, deltaX, deltaT, transitionLength, ctimestep, Tmp, L_curv,...
            endTime, seN, N_Reals, startRealization, cutOff_thickness, R_f, Rc, correctionLP_switch, repulsion_coeff, repulsion_expon, vdW_repulsion, res_limit, eq_thickness_EDL_vdW)
        cd ..
    case 'finiteSizedNonFlatFilms2D'
        mkdir(filmConfiguration)
        cd(filmConfiguration)
        %% raw conditions from the paper
%         R_film = [20 25 30 35 40 50 60 65 70 75 80 85 90 100 115 150 200 300 400 500 600 700 800 900 1000];   % radius of the film
        R_film = [100];
        R_f = R_film.*10^-6;              % in m
        h0_init = 2000e-9;                % initial film height in m
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
        h_critical_start = 0.627*kappa^(-2/7)*1.2      % this is when the film thinning velocity starts becoming nearly constant
        h_critical_end = 0.627*kappa^(-2/7)*0.8        % this is when the film thinning velocity ends becoming nearly constant
        repulsion_coeff = 1;
        repulsion_expon = 1;
        vdW_repulsion = 1;
        eq_thickness_EDL_vdW = 1;
        %% domain size and discretization parameters

        deltaX = 0.0025*ones(size(R_f));            % grid size (tested for grid independent results)
        L_flat = round(R_f/l_scale,4);              % length of the flat film
        N_flat = round(L_flat./deltaX);             % number of grid points in the same
%         ctimestep = 3;                              % exponent used in deciding deltaT = deltaX^c --> although c = 2.75 suffices, but a higher temp resolution enables more time stamps
        for i = 1:length(R_film)
            if R_film(i) < 50
                seN(i) = 20;                            % save every these many time steps
                ctimestep(i) = 3;
            elseif R_film(i) >= 50 && R_film(i) < 400
                seN(i) = 20;
                ctimestep(i) = 3.0;
            else
                seN(i) = 5;                            % save every these many time steps
                ctimestep(i) = 4;
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
        transitionLength = 1*lowerLimitOnL_curv;
        L_curv = 15*lowerLimitOnL_curv;                       % length of the curved portion of the film, for kappa > 1, one needs a smaller value of of L_curv
        endTime = 0.01;
        
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
            mkdir('disjPress_on')
            cd('disjPress_on')
        elseif isequal(disjPress_switch, 'off') 
            cutOff_thickness = 0.01;             % since the thinning rate in the absence of disj pres decreases asymptotically, a higher cut-off would save computational time
            mkdir('disjPress_off')
            cd('disjPress_off')
        end
        
        
        %% parent folder
        
        mk = strcat('h0_',num2str(h0_init*10^9),'nm','_Avw_',num2str(A_vw),'_ST_',num2str(gam),'_Rc_',num2str(Rc), '_disjPr_',disjPress_switch);
        mkdir(mk);
        cd(mk);
        copyfile('../../../*.m', '.')
        copyfile('../../../dataRadoev1984.xlsx', '.')
        
        %% plot the first film length
        [h x] = initialProfile(kappa,L_flat(1),L_curv,R_f, Rc, transitionLength,deltaX(1),filmConfiguration, correctionLP_switch);
        length(h)
        plot(x,h,'o')
        ylim([0.95 5])
        %% start simulations for different radii of films

        for i = 1:length(L_flat)
            str1{i} = strcat('Rf_', num2str(R_film(i)),'_mu_m');
            mkdir(str1{i});
            path_dest{i} = strcat('./', num2str(str1{i}));
            destinatn = path_dest{i};
            copyfile('*.m', destinatn)
            run_mainFiles{i} = strcat('./',num2str(str1{i}));
            cd(run_mainFiles{i})
            [t_rupt(i) drainageTime(i) drainageTime_right(i) drainageTime_right_rupt(i)  ...
            avg_cr_thinningRate_fit(i) h_cr_final(:,i) h_cr_final_FullFilmavg(:,i)] = main_finiteSizedFilms(filmConfiguration, disjPress_switch, ...
                        R_f(i), h0_init, A_vw, gam, Rc, visc, L_flat(i), N_flat(i), deltaX(i), deltaT(i), transitionLength, h_drain_start,...
                        h_drain_end, h_critical_start, h_critical_end, t_cr_dimensional, res_limit, ctimestep(i), Tmp, L_curv, endTime, seN(i),...
                        N_Reals, animationSkip, startRealization, hJoyeStart, hJoyeEnd, cutOff_thickness, correctionLP_switch, repulsion_coeff,...
                        repulsion_expon, vdW_repulsion, eq_thickness_EDL_vdW);

            
            cd ..
            fileToBeSaved = strcat('workspace_','Rf_',num2str(min(R_film)),'_to_',num2str(max(R_film)),'.mat');
            save(fileToBeSaved)
        end
        save('results_differentFilmSize.mat')
        
        %% post process finite sized films
        
        postProcessFiniteSizedFilms();
        
    case 'axisSymmetricFilm'
        mkdir(filmConfiguration)
        cd(filmConfiguration)
        switch ExpOrPar
            case 'experimental'

%                 R_film =  [40 50 60 65 70 75 80 85 90 100 115 150 200 300 400 500 600 700 800 900 1000 1500 2000 3000 4000];   % radius of the film
%                 R_film = [900 1000 1100 1200 1300 1400 1500 2000 3000 4000];
%                 R_film = [40 50 80 100 150 200 300 400 800 2000 4000];
                R_film = [4000];

                R_f = R_film.*10^-6;                            % in m
                h0_init = 2000e-9;                               % initial film height in m
                A_vw = 1.5e-20;                                % Hamaker constant
                gam = 0.0445;                                    % surface tension
                Rc = 1.8e-3;                                    % radius of capillary
                visc = 0.00089;                                 % viscosity
                t_cr_dimensional = [2 3 4 5];                   % time resolution (in sec) 
                res_limit_dimensional = [22e-6];  % spatial resolution (in m)

                %% derived quantities from above 

                l_scale = sqrt(Rc*h0_init/4);                   % length scale of the system - obtained from the O(1) scaling
                kappa = pi*h0_init^3*gam/A_vw/Rc;               % dimensionless curvature - the free parameter of the system
                t_scale = 3*visc*Rc^2/(2*gam*h0_init);            % time scale of the system

                %% drainage time start and end, critical thinning rates start and end, Joye's thinning rate start and end

                res_limit = res_limit_dimensional/l_scale;          % spatial resolution (scaled)
                h_drain_start = 100e-9/h0_init;                     % 100nm as mentioned in Wasan & Malhotra (1987)
                h_drain_end = [50e-9 40e-9 30e-9 25e-9]/h0_init;                        % 25 nm as mentioned in Wasan & Malhotra (1987)
%                 h_drain_end = [25e-9]/h0_init;                        % 25 nm as mentioned in Wasan & Malhotra (1987)
                hJoyeStart = 0.8;                                   % determine when to start measuring thinning rates
                hJoyeEnd = 0.627*kappa^(-2/7);                      % where to end
        %         h_critical_start = 0.627*kappa^(-2/7)*1.2;        % this is when the film thinning velocity starts becoming nearly constant
        %         h_critical_end = 0.627*kappa^(-2/7)*0.8;          % this is when the film thinning velocity ends becoming nearly constant
                h_critical_start = 0.627*kappa^(-2/7)*2;            % this is when the film thinning velocity starts becoming nearly constant
                h_critical_end = 0.627*kappa^(-2/7)*1.6;            % this is when the film thinning velocity ends becoming nearly constant

                %% domain size and discretization parameters

                deltaX = 0.05*ones(size(R_f));             % grid size (tested for grid independent results
                L_flat = round(R_f./l_scale,1);            % length of the flat film 
                N_flat = round(L_flat./deltaX);            % number of grid points in the same
                ctimestep = 2.25;                           % exponent used in deciding deltaT = deltaX^c --> although c = 2.75 suffices, but a higher temp resolution enables more time stamps
                seN = 100;                                 % save every seN time steps

                for i = 1:length(N_flat)
                    if (N_flat(i)) < 40
                        N_flat(i) = 40;
                        deltaX(i) = L_flat(i)/N_flat(i);
                    end        
                end

                %% simulation parameters

                Tmp = 0.0;                          % dimensionless noise strength (= 0, for deterministic)
                upperLimitOnL_curv = sqrt(2*Rc/h0_init) + sqrt(2*Rc/h0_init - L_flat.^2);
                lowerLimitOnL_curv = sqrt(h0_init*Rc)./l_scale;
%                 transitionLength = lowerLimitOnL_curv;
                transitionLength = 2;
                L_curv = 10*lowerLimitOnL_curv;                       % length of the curved portion of the film, for kappa > 1, one needs a smaller value of of L_curv
                endTime = 6000;

                N_Reals = 1;                        % number of realizations
                for i = 1:length(R_film)
                    if R_film(i) < 50
                        ctimestep(i) = 2.0;
                        seN(i) = 20.0;                            % save every these many time steps
                        animationSkip(i) = 10;                 % save animation every these many time steps
                    elseif R_film(i) >= 50 && R_film(i) <= 100
                        ctimestep(i) = 2.0;
                        seN(i) = 20;                            % save every these many time steps
                        animationSkip(i) = 10;                  % save animation every these many time steps      
                    elseif R_film(i) >= 900
                        ctimestep(i) = 2.0;
                        seN(i) = 20;                            % save every these many time steps
                        animationSkip(i) = 200; 
                    elseif R_film(i) == 4000 || R_film(i) == 50000 
                        ctimestep(i) = 2.0;
                        seN(i) = 20;                            % save every these many time steps
                        animationSkip(i) = 200; 
                    else
                        ctimestep(i) = 2.0;
                        seN(i) = 20;                            % save every these many time steps
                        animationSkip(i) = 200;                  % save animation every these many time steps
                    end
                    
                end
                deltaT = deltaX.^ctimestep;

                startRealization = 1;               % first realization

        %         for i = 1:length(L_flat)
        %             if L_flat(i) <= transitionLength
        %                 error('Length of the flat film is smaller than the transition region')
        %                 break;
        %             end
        %         end

                if isequal(disjPress_switch, 'on') 
                    cutOff_thickness = 1e-05;            % keep it lower to be able to probe even smaller thicknesses if it reaches
                    mkdir('disjPress_on')
                    cd('disjPress_on')
                    if isequal(repulsion_switch,'on')
        %                 mkdir('repulsion_on')
        %                 cd('repulsion_on')
                        if isequal(adhocRepulsion_switch, 'on')
                            mkdir('repulsion_adhoc_on')
                            cd('repulsion_adhoc_on')
                            c1 = 0;
                            c2 = 0;
                            c3 = 5*10^-9/h0_init;  
                            cutOff_thickness = h_drain_end;
                            eq_thickness_EDL_vdW = 0.0;

                            %% parent folder

                            mk = strcat('adj_CP_h0_',num2str(h0_init*10^9),'nm','_Avw_',num2str(A_vw),'_ST_',num2str(gam),'_Rc_',num2str(Rc), ...
                                        '_100_25nm','_c1_',num2str(c1),'_c2_',num2str(c2),'_hmin_',num2str(c3*h0_init*10^9),'nm');
                            mkdir(mk);
                            cd(mk);
                            copyfile('../../../../*.m', '.')
                            copyfile('../../../../dataRadoev1984.xlsx', '.')
                        elseif isequal(adhocRepulsion_switch, 'off')
                            mkdir('repulsion_edl_on')
                            cd('repulsion_edl_on')
                            c1 = 80;
                            c2 = 50;
            %                 c1 = 494.44;
            %                 c2 = 6e6;
                            c3 = 0;
                            eq_thickness_EDL_vdW = 0.065;           % you can take this value from the disj pressure isotherm.. I should rather call a function that calculates this for different choices of A_vdW, c1 and c2.
                            cutOff_thickness = h_drain_end;
                            %% parent folder

                            mk = strcat('adj_CP_h0_',num2str(h0_init*10^9),'nm','_Avw_',num2str(A_vw),'_ST_',num2str(gam),'_Rc_',num2str(Rc), ...
                                        '_100_25nm','_c1_',num2str(c1),'_c2_',num2str(c2));
                            mkdir(mk);
                            cd(mk);
                            copyfile('../../../../*.m', '.')
                            copyfile('../../../../dataRadoev1984.xlsx', '.')
                        end
                    elseif isequal(repulsion_switch,'off')
                        mkdir('repulsion_off')
                        cd('repulsion_off')
                        c1 = 0;
                        c2 = 0;
                        c3 = 0;
                        eq_thickness_EDL_vdW = 0;
                        mk = strcat('h0_',num2str(h0_init*10^9),'nm','_Avw_',num2str(A_vw),'_ST_',num2str(gam),'_Rc_',num2str(Rc), ...
                                        '_100_25nm','_c1_',num2str(c1),'_c2_',num2str(c2),'clr_One');
                        mkdir(mk);
                        cd(mk);
                        copyfile('../../../../*.m', '.');
                        copyfile('../../../../dataRadoev1984.xlsx', '.');
                    end         
                elseif isequal(disjPress_switch, 'off') 
                    cutOff_thickness = h_drain_end;             % since the thinning rate in the absence of disj pres decreases asymptotically, a higher cut-off would save computational time
                    mkdir('disjPress_off')
                    cd('disjPress_off')
                    mkdir('repulsion_off')
                    cd('repulsion_off')
                    c1 = 0;
                    c2 = 0;
                    c3 = 0;
                    eq_thickness_EDL_vdW = 0;
                    mk = strcat('Unadj_CP_h0',num2str(h0_init*10^9),'nm','_Avw_',num2str(A_vw),'_ST_',num2str(gam),'_Rc_',num2str(Rc), ...
                                '_100_25nm','_c1_',num2str(c1),'_c2_',num2str(c2),'_June13th_no_vdW');
                    mkdir(mk);
                    cd(mk);
                    copyfile('../../../../*.m', '.')
                    copyfile('../../../../dataRadoev1984.xlsx', '.')

                end
                
                 %% plot the first film length
                hfig = figure;
                hfig.Renderer = 'Painters';
                figureName = strcat('initialProfile_h0_', num2str(h0_init*1e9), 'nm');
                [h x] = initialProfile(kappa,L_flat(1),L_curv, R_f(1), Rc, transitionLength,deltaX(1),filmConfiguration, correctionLP_switch);
                area([0; x],[1; h])
                xlabel('$R_{film}$ [-]','Fontsize',14)
                ylabel('$h$ [-]','Fontsize',14)
                set(gca,'FontSize',14)
                ylim([0 2])



                 %% plot the first film length
                hfig = figure;
                hfig.Renderer = 'Painters';
                figureName = strcat('initialProfile_h0_', num2str(h0_init*1e9), 'nm');
                [h x] = initialProfile(kappa,L_flat(1),L_curv, R_f(1), Rc, transitionLength,deltaX(1),filmConfiguration, correctionLP_switch);
                area([0; x]*l_scale*10^6,[1; h]*h0_init*10^9)
                xlabel('$R_{film}$ ($\mu$m)','Fontsize',14)
                ylabel('$h$ (nm)','Fontsize',14)
                set(gca,'FontSize',14)
                ylim([0 2*h0_init*10^9])

                set(hfig,'Units','Inches');
                pos = get(hfig,'Position');
                set(hfig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
                print(hfig,figureName,'-dpdf','-r300');

                %% start simulations for different radii of films

              for i = 1:length(L_flat)
                    str1{i} = strcat('Rf_', num2str(R_film(i)),'_mu_m');
                    mkdir(str1{i});
                    path_dest{i} = strcat('./', num2str(str1{i}));
                    destinatn = path_dest{i};
                    copyfile('*.m', destinatn);
                    run_mainFiles{i} = strcat('./',num2str(str1{i}));
                    cd(run_mainFiles{i});
                    [t_rupt(i) drainageTime{i} drainageTime_right{i} drainageTime_right_rupt{i}, ...
                    avg_cr_thinningRate_fit(i) h_cr_final(:,i) h_cr_final_FullFilmavg(:,i) R_film_final(i) Press{i} vol_end(i) area_end(i) ...
                    ratio_vol(i) ratio_area(i) t_series_rel{i} h_fit_MTR{i} bb(:,i) h_c_tr(:,i)] = main_axisSymmetryFilm(filmConfiguration, disjPress_switch, correctionLP_switch, ...
                                R_f(i), h0_init, A_vw, c1, c2, c3, gam, Rc, visc, L_flat(i), N_flat(i), deltaX(i), deltaT(i), transitionLength, h_drain_start,...
                                h_drain_end, h_critical_start, h_critical_end, t_cr_dimensional, res_limit, ctimestep(i), Tmp, L_curv, endTime, seN(i),...
                                N_Reals, animationSkip(i), startRealization, hJoyeStart, hJoyeEnd, cutOff_thickness, eq_thickness_EDL_vdW);


                    cd ..
                    fileToBeSaved = strcat('workspace_','Rf_',num2str(min(R_film)),'_to_',num2str(max(R_film)),'.mat');
                    save(fileToBeSaved)
                end
%                 save(strcat('results_differentFilmSize_diffThinnRates',num2str(R_f),'.mat'))
                close all;
                save(strcat('results_differentFilmSize_diffThinnRates_04','.mat'))

                %% post process finite sized films

                postProcessFiniteSizedFilms();
            case 'parametric'
                mkdir(ExpOrPar)
                cd(ExpOrPar)
                kappa = [10^-5; 10^-4; 10^-3; 10^-2; 10^-1; 1; 10; 10^2; 10^3];
                Tmp = 0.0;
                L_flat = 5.7; 
                R_f = 1;
                Rc = 1;
                res_limit = 22e-6;
                transitionLength = 1./sqrt(2*kappa);
                L_curv = 10*transitionLength;

                if Tmp == 0
                    N_Reals = 1;
                else
                    N_Reals = 2;   %% please adjust the number of realizations based on how many you want to sample
                end

                if isequal(disjPress_switch, 'on') 
                    cutOff_thickness = 1e-05;            % keep it lower to be able to probe even smaller thicknesses if it reaches 
                    mkdir('disjPress_on')
                    cd('disjPress_on')
                    if isequal(repulsion_switch,'on')
                        if isequal(adhocRepulsion_switch, 'on')
                            mkdir('repulsion_adhoc_on')
                            cd('repulsion_adhoc_on')
                            c1 = 0;
                            c2 = 0;
                            c3 = 5*10^-9/h0_init;  
                            cutOff_thickness = h_drain_end;
                            eq_thickness_EDL_vdW = 0.0;

                            %% parent folder
                            copyfile('../../../../../*.m', '.')
                            copyfile('../../../../../dataRadoev1984.xlsx', '.')
                        elseif isequal(adhocRepulsion_switch, 'off')
                            mkdir('repulsion_edl_on')
                            cd('repulsion_edl_on')
                            c1 = 80;
                            c2 = 50;
            %                 c1 = 494.44;
            %                 c2 = 6e6;
                            c3 = 0;
                            eq_thickness_EDL_vdW = 0.065;           % you can take this value from the disj pressure isotherm.. I should rather call a function that calculates this for different choices of A_vdW, c1 and c2.
                            cutOff_thickness = h_drain_end;
                            %% parent folder
                            copyfile('../../../../*.m', '.')
                            copyfile('../../../../dataRadoev1984.xlsx', '.')
                        end
                    elseif isequal(repulsion_switch,'off')
                        mkdir('repulsion_off')
                        cd('repulsion_off')
                        c1 = 0;
                        c2 = 0;
                        c3 = 0;
                        eq_thickness_EDL_vdW = 0;
                        
                    end         
                elseif isequal(disjPress_switch, 'off') 
                    cutOff_thickness = 0.05;             % since the thinning rate in the absence of disj pres decreases asymptotically, a higher cut-off would save computational time
                    mkdir('disjPress_off')
                    cd('disjPress_off')
                end

                deltaX = 0.05;
                N = (L_flat+L_curv)./deltaX;
                ctimestep = 2.5;
                deltaT = deltaX^ctimestep;
                endTime = 20000;
                seN = 50;
                startRealization = 1;

                [h x] = initialProfile(kappa,L_flat(1),L_curv(1), R_f(1), Rc, transitionLength,deltaX(1),filmConfiguration, correctionLP_switch);

                %% here we kick-start the simulations

                str1 = strcat('theta_', num2str(Tmp));
                mkdir(str1);
                cd(str1);
                copyfile('../../../../../*.m', '.')
                for i = 1:length(kappa)
                    str2{i} = strcat('kappa_', num2str(kappa(i)));
                    mkdir(str2{i});
                    path_dest = strcat('./', num2str(str2{i}));
                    destinatn = path_dest;
                    copyfile('../*.m', destinatn)
                    [t_rupt(i) x_rupt(i)] = main_axisSymmetryFilm_parametric(filmConfiguration,disjPress_switch, kappa, L_flat, deltaX, deltaT, transitionLength, ctimestep, Tmp, L_curv,...
                    endTime, seN, N_Reals, startRealization, cutOff_thickness, correctionLP_switch, c1, c2, c3, R_f, Rc, res_limit,...
                                                             eq_thickness_EDL_vdW)
                end
                            

                
        end
        
end

save('workspaceInputParameters.mat')

% cd ../../../../

toc

%%
%%%%%%%%%%%%%%%%%%%%%% end this file here  %%%%%%%%%%%%%%%%%%%%%%%%%%%




