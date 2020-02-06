function [h_right_avg_j x_dimple_loc_right h_right x_right vol area_film] = spatialResolution_filmThickness(x,c,t,res_limit, deltaX);
% function [h_right_avg_j h_left_avg_j x_dimple_loc_right h_right x_right] = spatialResolution_filmThickness(x,c,t,res_limit, deltaX);

%% calculate region around the right dimple
%     x = x';
    x_ri = x(x>0);
    for i = 1:length(t);
        h_r = c(:,i);
        h_ri(:,i) = h_r(x>0);
%         sizehri = size(h_ri)
%         [idx_rightMin(i)] = max(find(round(h_ri(:,i),3) == min(round(h_ri(:,i),3))));
        [idx_rightMin(i)] = max(find(h_ri(:,i) == min(h_ri(:,i))));
        if mod(round(res_limit/deltaX),2) ~= 0
            if ( x_ri(idx_rightMin(i)) - x_ri(round((res_limit/2)/deltaX)) ) > 0              
                h_right{i} = h_ri((x_ri > (x_ri(idx_rightMin(i) - round((res_limit/2)/deltaX))) & x_ri <= (x_ri(idx_rightMin(i) + round((res_limit/2)/deltaX)))),i);
                x_right{i} = x_ri( x_ri > (x_ri(idx_rightMin(i) - round((res_limit/2)/deltaX))) & x_ri <= (x_ri(idx_rightMin(i) + round((res_limit/2)/deltaX))))';
            else
                h_right{i} = h_ri((x_ri > (x_ri(idx_rightMin(i))) & x_ri <= (x_ri(idx_rightMin(i) + round(res_limit/deltaX) + 1))),i);
                x_right{i} = x_ri( x_ri > (x_ri(idx_rightMin(i))) & x_ri <= (x_ri(idx_rightMin(i) + round(res_limit/deltaX) + 1)))';
            end
        else
            if ( x_ri(idx_rightMin(i)) - x_ri(round((res_limit/2)/deltaX)) ) > 0              
                h_right{i} = h_ri((x_ri > (x_ri(idx_rightMin(i) - round((res_limit/2)/deltaX))) & x_ri <= (x_ri(idx_rightMin(i) + round((res_limit/2)/deltaX)))),i);
                x_right{i} = x_ri( x_ri > (x_ri(idx_rightMin(i) - round((res_limit/2)/deltaX))) & x_ri <= (x_ri(idx_rightMin(i) + round((res_limit/2)/deltaX))))';
            else
                h_right{i} = h_ri((x_ri > (x_ri(idx_rightMin(i))) & x_ri <= (x_ri(idx_rightMin(i) + round(res_limit/deltaX) ))),i);
                x_right{i} = x_ri( x_ri > (x_ri(idx_rightMin(i))) & x_ri <= (x_ri(idx_rightMin(i) + round(res_limit/deltaX) )))';
            end
        end
            
        h_right_avg(i) = mean(h_right{i});
        x_dimple_loc_right(i) = x_ri(idx_rightMin(i));
%         x_ri(1:idx_rightMin(i))
%         h_ri(1:idx_rightMin(i))
%         aa0 = size(h_ri(1:idx_rightMin(i)),i)
%         aa = size(h_ri(1:idx_rightMin(i),i).*x_ri(1:idx_rightMin(i)))
%         bb = size(x_ri(1:idx_rightMin(i)))
        if idx_rightMin(i) ~= 1
            area_film(i) = trapz(x_ri(1:idx_rightMin(i)), h_ri(1:idx_rightMin(i),i));
            vol(i) = 2*pi*trapz(x_ri(1:idx_rightMin(i)), h_ri(1:idx_rightMin(i),i).*x_ri(1:idx_rightMin(i)));
        else
            area_film(i) = 1; %area_film(i-1);
            vol(i) = 1; %vol(i-1);
        end
    end
    h_right_avg_j = h_right_avg(:);
    h_right_avg_j(h_right_avg_j == 0) = [];
    
%     idx_rightMin
%     x_ri(1:idx_rightMin)
%     h_ri(1:idx_rightMin)
%     vol = trapz(h_ri(1:idx_rightMin).*x_ri(1:idx_rightMin), x_ri(1:idx_rightMin));
    
    %% dummy
%     x_le = x(x<0);
%     x_le = flipud(x_le);
%     for i = 1:length(t)-1;
%         h_l = c(:,i);
%         h_le(:,i) = h_l(x<0);
%         h_le(:,i) = flipud(h_le(:,i));
%         [val_leftMin(i) idx_leftMin(i)] = min(h_le(:,i));
%         if mod(round(res_limit/deltaX),2) ~= 0
%             if ( x_le(idx_leftMin(i)) - x_le(round((res_limit/2)/deltaX)) ) < 0              
%                 h_left{i} = h_le( x_le < (x_le(idx_leftMin(i) - round((res_limit/2)/deltaX))) & x_le >= (x_le(idx_leftMin(i) + round((res_limit/2)/deltaX))),i);
%                 x_left{i} = x_le( x_le < (x_le(idx_leftMin(i) - round((res_limit/2)/deltaX))) & x_le >= (x_le(idx_leftMin(i) + round((res_limit/2)/deltaX))))';
%             else
%                 h_left{i} = h_le((x_le < (x_le(idx_leftMin(i))) & x_le >= (x_le(idx_leftMin(i) + round(res_limit/deltaX) + 1))),i);
%                 x_left{i} = x_le( x_le < (x_le(idx_leftMin(i))) & x_le >= (x_le(idx_leftMin(i) + round(res_limit/deltaX) + 1)))';
%             end
%         else
%             if ( x_le(idx_leftMin(i)) - x_le(round((res_limit/2)/deltaX)) ) < 0              
%                 h_left{i} = h_le( x_le < (x_le(idx_leftMin(i) - round((res_limit/2)/deltaX))) & x_le >= (x_le(idx_leftMin(i) + round((res_limit/2)/deltaX))),i);
%                 x_left{i} = x_le( x_le < (x_le(idx_leftMin(i) - round((res_limit/2)/deltaX))) & x_le >= (x_le(idx_leftMin(i) + round((res_limit/2)/deltaX))))';
%             else
%                 h_left{i} = h_le((x_le < (x_le(idx_leftMin(i))) & x_le >= (x_le(idx_leftMin(i) + round(res_limit/deltaX) ))),i);
%                 x_left{i} = x_le( x_le < (x_le(idx_leftMin(i))) & x_le >= (x_le(idx_leftMin(i) + round(res_limit/deltaX) )))';
%             end
%         end
%             
%         h_left_avg(i) = mean(h_left{i});
%         x_dimple_loc_left(i) = x_le(idx_leftMin(i));
%     end
%     h_left_avg_j = h_left_avg(:);
%     h_left_avg_j(h_left_avg_j == 0) = [];
    
end