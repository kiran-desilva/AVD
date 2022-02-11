

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
%%		-	semispan
%%
%%				Material Parameters
%%--------------------------------------------------------------
%%		All in metric
%%
%%		-	Young's modulus (E when a variable)
%%
%% 				Design Parameters
%%--------------------------------------------------------------
%%		-	stringer pitch
%%		-	stringer thickness
%%		-	stringer web height
%%		-	stringer flange to web ratio
%%
%%					Other
%%--------------------------------------------------------------
%%		- spanwise bending moment distribution
%%
%%				Return Value
%%		- The total volume used up by stringers + skin
%%		  this value is proportional to weight and thus
%%		  should be minimised

function total_volume = rib_stringer_func(geometry, material, design_params, bending_moment_dist)
	spanwise_station = 0; % Start at root
	total_volume = 0;

	while spawise_station < geometry.semispan
		bending_moment = double(bending_moment_dist(spanwise_station));

		box_width = double(geometry.box_width_func(spanwise_station));
		web_height = double(geometry.web_height_func(spanwise_station));

		stringer.flange_width = design_params.stringer_web_height*design_params.flange_to_web_ratio;
		stringer.cross_sec_area = ...; % TODO: As appropriate for the selected stringer type

		K = 0.9*4; % From buckling of SS plate plot, taken for a/b -> infinity
		comp_load_per_length = bending_moment/(box_width*web_height);
		panel_thickness = (comp_load_per_length*design_params.stringer_pitch^2/(K*material.E))^(1/3);
		sigma_0 = comp_load_per_length/panel_thickness; % Critical Buckling Stress -> need to adjust this
		% for stringers using catchpole diagram

		eff_thickness = panel_thickness + stringer.cross_sec_area/design_params.stringer_pitch; 
		number_of_panels = ...;
		panel_eff_area = eff_thickness*design_params.stringer_pitch*number_of_panels;

	end
end