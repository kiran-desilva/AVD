function [updated_load_dist] = add_point_load(load_dist,x_dist,load,cg)
    %find closest station
    [~,idx] = min(abs(x_dist - cg));
    %update weight
    updated_load_dist = load_dist;
    updated_load_dist(idx) = updated_load_dist(idx) + load;
