function [struct] = wing_geometry_calc_struct(struct,n)
    [struct.Croot,struct.Ctip,struct.Cmac,struct.b,struct.Ybar,struct.Xac] = wing_geometry_calc(struct.s,struct.Ar,struct.lambda,struct.sweepLE,n);
end