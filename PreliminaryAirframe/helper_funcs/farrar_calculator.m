function F = farrar_calculator(A_s_over_bt, t_s_over_t)
%     if false
%         close all
%         %clc
%         
%         digitise_farrar;
%         hold on;
%         scatter(A_s_over_bt, t_s_over_t, 300, 'g', 'x');
%     end
    persistent fd;
    if isempty(fd)
        load farrar_data.mat farrar_data;
        fd = farrar_data;
    end
    
    d = fd.transform([A_s_over_bt, t_s_over_t]);
    [theta, ~] = cart2pol(d(:, 1), d(:, 2));
    distances = zeros(1, numel(fd.fit));
    for i = 1:numel(fd.fit)
        XY = fd.fit{i}(theta);
        distances(i) = sqrt((XY(:, 1) - A_s_over_bt).^2 + (XY(:, 2) - t_s_over_t).^2);
        % scatter(XY(:, 1), XY(:, 2));
    end
    [best_fits, best_fit_idx] = sort(distances);
    distance_between_best = abs(best_fits(1) + best_fits(2));
    F = (best_fits(1)*fd.f_vals(best_fit_idx(1)) + best_fits(2)*fd.f_vals(best_fit_idx(2)))/distance_between_best;
    %hold off;
end

function F = weigh_by_distance(curve_distances, f_vals)
    [sorted_distance, sorted_idx] = sort(curve_distances);
    total_distance = sum(abs(sorted_distance));
    F = sum(sorted_distance.*f_vals(sorted_idx))/total_distance;
end