function [struct] = wing_geometry_calc_struct(struct,n)
    [struct.Croot,struct.Ctip,struct.Cmac,struct.b] = wing_geometry_calc(struct.s,struct.Ar,struct.lambda,struct.sweep_25,n);
end