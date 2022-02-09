function [load,cg] = point_load_from_dist(x,loaddist)
    load = sum(loaddist);
    cg = sum(loaddist.*x)/load;
