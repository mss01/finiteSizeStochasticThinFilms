function [h_right_avg_j h_min h_centre_j v_thin_rim v_thin_centre avg_cr_thinningRate_fit h_cr_final ...
    h_cr_final_FullFilmavg drainageTime drainageTime_right drainageTime_right_rupt t_rupt h_max_dimp_r...
    beginDrainageTime_right endDrainageTime_right x_dimple_loc_right vol_end area_end h_min_fit_AtDimple] = loopingOverRealizations(x, h0_init, L_flat, transitionLength, tt, x_centre, ...
                                        locDimple_right, res_limit, cr_thickness, h_drain_start, h_drain_end, ...
                                        t_cr, deltaX, deltaT, h_critical_start, h_critical_end);
                                    

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaulttextInterpreter','latex');
                              
% newFolder = 'realization1';
% str1 = strcat('.\',newFolder, '\*.txt');
% textFiles = dir([str1]);
% a   = cellfun(@num2str, struct2cell(textFiles), 'UniformOutput', false);
% Out = sortrows(a.',6);
load('realization1/hData.mat');
b = length(t_store);
t = t_store;                    % corresponds to each time step, including the initial time step                    

%% read the height profiles for every time steps
% for k = 1:1:b   % 1:s:(b-last)
%     str2 = strcat('.\', newFolder, '\', Out{k,1});
%     c(:,k)  = dlmread(str2);
% end
% c(:,all(c==0)) = [];
% c = [c(:,(1:end))];
q = max(size(t));
for i = 1:q-1
    Y = h_store(:,i);
    h_centre(i) = Y(x_centre);
%     h_edges_left(i) = Y(locDimple_left);
    h_edges_right(i) = Y(locDimple_right);
%     h_avg(i) = mean(Y(locDimple_left:x_centre));
    h_avg(i) = mean(Y(x_centre:locDimple_right));
    t_centre(i) = t(i+1);
%     Y1 = h_store((x > (-L_flat + transitionLength)& x <= 0),i);
    Y2 = h_store((x >= 0 & x < (L_flat)),i);
%     Y2 = h_store((x >= 0 & x < (L_flat - transitionLength)),i);
%     h_max_dimp_l(i) = max(Y1);
    [h_max_dimp_r(i) idx_h_max(i)] = max(Y2(1:end-1));
    x_max_dimp_r(i) = x(idx_h_max(i));
    if idx_h_max(i) == 1;
        Pres_h_max(i) = -1./(x(idx_h_max(i)+1).*deltaX^2).*((x(idx_h_max(i)+1) + x(idx_h_max(i)+2))./2.*Y2((idx_h_max(i)+2)) - ((x(idx_h_max(i)+1) + x(idx_h_max(i)+2))./2 + (x(idx_h_max(i)+1) + x(idx_h_max(i)))./2).*Y2(idx_h_max(i)+1) + ...
                    (x(idx_h_max(i)+1) + x(idx_h_max(i)))./2.*Y2(idx_h_max(i)))./deltaX^2;
    else
        Pres_h_max(i) = -1./(x(idx_h_max(i)).*deltaX^2).*((x(idx_h_max(i)) + x(idx_h_max(i)+1))./2.*Y2(idx_h_max(i)+1) - ((x(idx_h_max(i)) + x(idx_h_max(i)+1))./2 + (x(idx_h_max(i)) + x(idx_h_max(i)-1))./2).*Y2(idx_h_max(i)) + ...
                    (x(idx_h_max(i)) + x(idx_h_max(i)-1))./2.*Y2(idx_h_max(i)-1))./deltaX^2;
    end
    h_max(i) = max(Y);
    h_min(i) = min(Y);
    [h_min_dimp_r(i) idx_h_min(i)] = min(Y(1:end-1));
%     idx_neighbours(:,i) = [idx_h_min(i)-3 idx_h_min(i)-2 idx_h_min(i)-1 idx_h_min(i) idx_h_min(i)+1 idx_h_min(i)+2 idx_h_min(i)+3];
%     idx_neighbours(:,i) = [idx_h_min(i)-1:1:idx_h_min(i)+1];
%     x_dimple_neighbours(:,i) = x(idx_neighbours(:,i));
%     h_min_dimple_neighbours(:,i) = Y(idx_neighbours(:,i));
%     par_fit_dimple(:,i) = polyfit(x_dimple_neighbours, h_min_dimple_neighbours, 2);
%     x_dimple_fit(:,i) = [x_dimple_neighbours(1,i):0.005:x_dimple_neighbours(end,i)];
%     h_min_dimple_fit(:,i) = polyval(par_fit_dimple(:,i), x_dimple_fit(:,i));
%     h_min_fit_AtDimple(i) = min(h_min_dimple_fit(:,i));
    h_min_fit_AtDimple(i) = h_min(i);
    del(i,1) = min(Y);
end

t_rupt = t(end);
h_centre_j = h_centre(:);
h_centre_j(h_centre_j == 0) = [];
% h_edges_left_j =  h_edges_left(:);
% h_edges_left_j(h_edges_left_j == 0) = [];
h_edges_right_j =  h_edges_right(:);
h_edges_right_j(h_edges_right_j == 0) = [];

%% lets start first with comparison with Joye

v_thin_rim = 1./t_rupt;
v_thin_centre = (h_centre_j(end) - h_centre_j(1))./t_rupt;


%% calculate region around the dimple
% for j = 1:length(res_limit)
%     [h_right_avg_j(:,j) h_left_avg_j(:,j) x_dimple_loc_right(:,j)] = spatialResolution_filmThickness(x,c,t,res_limit(j));
% end
for i = 1:length(res_limit)
    [h_right_avg_j(:,i) x_dimple_loc_right(:,i) h_right(:,i) x_right(:,i) vol(:,i) area_film(:,i)] = spatialResolution_filmThickness(x,h_store,t,res_limit(i), deltaX);
end
% [h_right_avg_j h_left_avg_j x_dimple_loc_right h_right x_right] = spatialResolution_filmThickness(x,h_store,t,res_limit, deltaX);

vol_end = vol(end);
area_end = area_film(end);

for i = 1:length(h_max_dimp_r)
        curvOverall(i) = (h_max_dimp_r(i) - h_min(i))./(x_max_dimp_r(i) - x_dimple_loc_right(i)).^2;
end

hfig_hMinfit = figure;
hfig_hMinfit.Renderer = 'Painters';

plot(t(2:end), h_min_fit_AtDimple)
hold on
plot(t(2:end), h_min)
xlabel('t [-]','Fontsize',14)
ylabel('$h_{min}$ [-]','Fontsize',14)
set(gca,'FontSize',16)

set(hfig_hMinfit,'Units','Inches');
pos = get(hfig_hMinfit,'Position');
set(hfig_hMinfit,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hfig_hMinfit,'hMinFitVsT','-dpdf','-r300')
savefig('hMinFitVsT')

hfig_vol = figure;
hfig_vol.Renderer = 'Painters';

plot(t, vol, 'o')
xlabel('t [-]','Fontsize',14)
ylabel('vol [-]','Fontsize',14)
set(gca,'FontSize',16)

set(hfig_vol,'Units','Inches');
pos = get(hfig_vol,'Position');
set(hfig_vol,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hfig_vol,'volVsT','-dpdf','-r300')
savefig('volVsT')

hfig_area = figure;
hfig_area.Renderer = 'Painters';

plot(t, area_film, 'o')
xlabel('t [-]','Fontsize',14)
ylabel('area [-]','Fontsize',14)
set(gca,'FontSize',16)

set(hfig_area,'Units','Inches');
pos = get(hfig_area,'Position');
set(hfig_area,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hfig_area,'areaVsT','-dpdf','-r300')
savefig('areaVsT')

close all;
save('hProbeAvgLite.mat','Pres_h_max','x_min','x_max_dimp_r','h_max_dimp_r')
save('hProbeAvg.mat')



%% calculate thinning rates based on the min film height

[avg_cr_thinningRate_fit intercept_cr_thinningRate avg_cr_thinningRate_mean] = critical_dhdt(t, del, h_critical_start, h_critical_end);


%% to calculate critical thickness based on when the average film height goes below the cross over thickness
h_avg_j = h_avg(:);
h_avg_j(h_avg_j == 0) = [];

for j = 1:length(t_cr)
    [h_cr(:,j) h_cr_right(:,j)] = criticalThicknesses_tempRes(deltaT, t, h_avg_j, del, cr_thickness, h_right_avg_j, t_cr(j));
%     [h_cr(:,j) h_cr_left(:,j) h_cr_right(:,j)] = criticalThicknesses_tempRes(deltaT, t, h_avg_j, del, cr_thickness, h_right_avg_j, h_left_avg_j, t_cr(j));
%     h_cr((h_cr == 0))             = [];
%     h_cr_left((h_cr_left(:,j) == 0))   = [];
%     h_cr_right((h_cr_right(:,j) == 0)) = [];
%     if isempty(h_cr)
%         clear h_cr
%         h_cr_final(j) = (h_cr_left(j) + h_cr_right(j))/2;
%     elseif isempty(h_cr_left)
%         clear h_cr_left
%         clear h_cr_right
%         h_cr_final(j) = h_cr(j);
%     end
    h_cr_final(j) = h_cr_right(j);
    h_cr_final_FullFilmavg(j) = h_cr(j);
end

%% to calculate the drainage time when the average height goes below 200 nm and reaches 50 nm
for ith = 1:length(h_drain_end)
    for j = 1:length(res_limit)
        [drainageTime(ith,j) drainageTime_right(ith,j) drainageTime_right_rupt(ith,j) beginDrainageTime_right(ith,j) endDrainageTime_right(ith,j)] = ...
                                                            drainageTimes(del,t,h_avg_j, h_right_avg_j(:,j), h_drain_start, h_drain_end(ith), cr_thickness);
%     [drainageTime(j) drainageTime_right(j) drainageTime_left(j) drainageTime_right_rupt(j) drainageTime_left_rupt(j)] = ...
%         drainageTimes(del,t,h_avg_j, h_right_avg_j(:,j), h_left_avg_j(:,j), h_drain_start, h_drain_end, cr_thickness);
    end
end

save('pressure_atDimple.mat', 'Pres_h_max')


hfig1 = figure;
plot(t(2:end), h_avg);
hold on
plot(t, h_right_avg_j)
yli1 = yline(h_drain_start,'-.k','h(t_d = 0)','LineWidth',1.5);
yli.LabelHorizontalAlignment = 'center';
yli2 = yline(h_drain_end(1),':k','h(t_d)','LineWidth',1.5);
yli2.LabelHorizontalAlignment = 'center';
yli3 = yline(h_drain_end(2),':k','h(t_d)','LineWidth',1.5);
yli3.LabelHorizontalAlignment = 'center';
yli4 = yline(h_drain_end(3),':k','h(t_d)','LineWidth',1.5);
yli4.LabelHorizontalAlignment = 'center';
yli5 = yline(h_drain_end(4),':k','h(t_d)','LineWidth',1.5);
yli5.LabelHorizontalAlignment = 'center';
xlabel('t [-]','Fontsize',14)
ylabel('$\langle h \rangle$ [-]','Fontsize',14)
legend('$\langle h(R) \rangle$','$\langle h(probe) \rangle = 0 \mu m$','$\langle h(probe) \rangle = 11 \mu m$','$\langle h(probe) \rangle = 22 \mu m$','$\langle h(probe) \rangle = 44 \mu m$')
set(gca,'FontSize',14)

set(hfig1,'Units','Inches');
pos = get(hfig1,'Position');
set(hfig1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hfig1,'avgThicknessVsT','-dpdf','-r300')
savefig('avgThicknessVsT')

hfig2 = figure;
hfig2.Renderer = 'Painters';
plot(t(2:end), Pres_h_max);
xlabel('t [-]','Fontsize',14)
ylabel('P [-]','Fontsize',14)
set(gca,'FontSize',14)

set(hfig2,'Units','Inches');
pos = get(hfig2,'Position');
set(hfig2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hfig2,'Pres_hmax_VsT','-dpdf','-r300')
savefig('Pres_hmax_VsT')

% clear c
clear dhmindt
% clear del
clear h_cr_all_left
clear h_cr_all_right
clear('t_drain', 't_drain_left', 't_drain_right')



end