function [h_right_avg_j h_left_avg_j h_avg h_min h_centre_j v_thin_rim v_thin_centre avg_cr_thinningRate_fit h_cr_final ...
    h_cr_final_FullFilmavg drainageTime drainageTime_left drainageTime_right drainageTime_right_rupt drainageTime_left_rupt t_rupt] = ...
                                        loopingOverRealizations(x, tt, x_centre, locDimple_left, ...
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
load('realization1\hData.mat');
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
    h_edges_left(i) = Y(locDimple_left);
    h_edges_right(i) = Y(locDimple_right);
    h_avg(i) = mean(Y(locDimple_left:x_centre));
    h_avg_right(i) = mean(Y(x_centre:locDimple_right));
    t_centre(i) = t(i+1);
    h_max(i) = max(Y);
    h_min(i) = min(Y);
    del(i,1) = min(Y);
end

t_rupt = t(end);
h_centre_j = h_centre(:);
h_centre_j(h_centre_j == 0) = [];
h_edges_left_j =  h_edges_left(:);
h_edges_left_j(h_edges_left_j == 0) = [];
h_edges_right_j =  h_edges_right(:);
h_edges_right_j(h_edges_right_j == 0) = [];

%% lets start first with comparison with Joye

v_thin_rim = 1./t_rupt;
v_thin_centre = (h_centre_j(end) - h_centre_j(1))./t_rupt;


%% calculate region around the dimple
% for j = 1:length(res_limit)
%     [h_right_avg_j(:,j) h_left_avg_j(:,j) x_dimple_loc_right(:,j)] = spatialResolution_filmThickness(x,c,t,res_limit(j));
% end

[h_right_avg_j h_left_avg_j x_dimple_loc_right h_right x_right] = spatialResolution_filmThickness(x,h_store,t,res_limit, deltaX);

save('hProbeAvg.mat')

%% calculate thinning rates based on the min film height

[avg_cr_thinningRate_fit intercept_cr_thinningRate avg_cr_thinningRate_mean] = critical_dhdt(t, del, h_critical_start, h_critical_end);


%% to calculate critical thickness based on when the average film height goes below the cross over thickness
h_avg_j = h_avg(:);
h_avg_j(h_avg_j == 0) = [];

for j = 1:length(t_cr)
    [h_cr(:,j) h_cr_left(:,j) h_cr_right(:,j)] = criticalThicknesses_tempRes(deltaT, t, h_avg_j, del, cr_thickness, h_right_avg_j, h_left_avg_j, t_cr(j));
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
    h_cr_final(j) = 0.5*(h_cr_left(j) + h_cr_right(j));
    h_cr_final_FullFilmavg(j) = h_cr(j);
end

%% to calculate the drainage time when the average height goes below 200 nm and reaches 50 nm
for j = 1:length(res_limit)
    [drainageTime(j) drainageTime_right(j) drainageTime_left(j) drainageTime_right_rupt(j) drainageTime_left_rupt(j)] = ...
        drainageTimes(del,t,h_avg_j, h_right_avg_j(:,j), h_left_avg_j(:,j), h_drain_start, h_drain_end, cr_thickness);
end

hfig1 = figure;
plot(t(2:end), h_avg);
hold on
plot(t(2:end), h_right_avg_j)
xlabel('time','Fontsize',14)
ylabel('$<h>$','Fontsize',14)
legend('$<h(R)>$','$<h(probe)>$')
set(gca,'FontSize',14)

set(hfig1,'Units','Inches');
pos = get(hfig1,'Position');
set(hfig1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(hfig1,'fullAndProbeaverageThicknessVsTime','-dpdf','-r300')

% clear c
clear dhmindt
% clear del
clear h_cr_all_left
clear h_cr_all_right
clear('t_drain', 't_drain_left', 't_drain_right')


end