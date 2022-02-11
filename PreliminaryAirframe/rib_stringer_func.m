

%% Places the ribs for a given geometry and free parameters
%% and calculates the total area, which can be used as an
%% optimisation metric to be minimised (for smallest weight
%% design)
%%
%%					Geometry
%%-------------------------------------------------------------
%%		-	chord function
%%		-	sweep angle
%%		-	spar placement as
%%			percentage of chord
%%		-	airfoil profile as
%%			a linear interp (NOT SURE IF THIS IS NEEDED
%%			it can maybe be replaced with just height/length
%%		 	of the wing box symfunc)
%%
%% 				Free Parameters
%%--------------------------------------------------------------
%%		-	number of panels at root (basically no. of stringers at root plus 1)
%%		-	stringer thickness
%%		-	stringer web height
%%		-	flange width of stringer
%%
%%					Other
%%--------------------------------------------------------------
%%		- spanwise bending moment distribution

function total_area = rib_stringer_func()

end