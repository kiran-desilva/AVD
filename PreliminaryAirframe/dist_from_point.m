function [dist] = dist_from_point(x_disc,load,cg)
    dist = zeros(size(x_disc));
    %find closest station
    [~,idx] = min(abs(x_disc - cg));
    %update weight
    dist(idx) = load;