function [load_dist] = add_point_load(load_dist,x_dist,load,cg)
    %find closest station
    [~,idx] = min(abs(x_dist - cg));
    %update weight
    load_dist(idx) = load_dist(idx) + load;
