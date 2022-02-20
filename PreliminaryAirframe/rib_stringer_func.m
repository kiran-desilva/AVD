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
%%		- 	Poissons ratio
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

function output = rib_stringer_func(geometry, material, design_params, bending_moment_dist, doPlot)
	spanwise_station = 0; % Start at root

	x_leading_edge = @(y) -tand(geometry.sweep_deg)*y + geometry.c(0)/2;
	x_trailing_edge = @(y) x_leading_edge(y) - geometry.c(y);
	x_at_some_percent_chord = @(y, percent_chord) x_leading_edge(y) - geometry.c(y)*percent_chord;
	x_front_spar = @(y) x_at_some_percent_chord(y, geometry.spar.front_x_c);
	x_rear_spar = @(y) x_at_some_percent_chord(y, geometry.spar.rear_x_c);

	options = optimoptions('fsolve','Display','none');

	starting_no_of_stringers = floor(geometry.box_width_func(0)/design_params.stringer_pitch);

	% Place the stringers s.t. the middle stringer connects middle of wingbox
	% at root to middle of wingbox at tip
	x0 = (x_front_spar(0) + x_rear_spar(0))/2;
	x1 = (x_front_spar(geometry.semispan) + x_rear_spar(geometry.semispan))/2;
	slope = (x1 - x0)/(geometry.semispan);

	stringer_x_space = linspace(0, 1, starting_no_of_stringers + 2);
	stringer_x_space = stringer_x_space(2:end-1);
	stringer_func = @(y, root_intercept_percent) slope*y + x_front_spar(0) - geometry.box_width_func(0)*root_intercept_percent;

	placed_stringer_func = @(y) stringer_func(y, stringer_x_space');
    intercepts = min([fsolve(@(x) x_front_spar(x) - placed_stringer_func(x), zeros(starting_no_of_stringers, 1), options),...
        fsolve(@(x) x_rear_spar(x) - placed_stringer_func(x), zeros(starting_no_of_stringers, 1), options)], [], 2);

	stringer_length_till_y = @(y, stringer_idx) vecnorm([y.*ones(size(stringer_idx)); stringer_func(y, stringer_x_space(stringer_idx))]...
                                                - [zeros(size(stringer_idx)); stringer_func(0, stringer_x_space(stringer_idx))]);
    percentage_of_stringer_at_y = @(y, stringer_idx) stringer_length_till_y(y, stringer_idx)./stringer_length_till_y(intercepts(stringer_idx)', stringer_idx);

	if doPlot
		y_space = linspace(0, geometry.semispan, 1000);
		figure;
		hold on;
		plot([0, 0], [x_leading_edge(0), x_trailing_edge(0)], 'k');
		plot([geometry.semispan, geometry.semispan], [x_leading_edge(geometry.semispan), x_trailing_edge(geometry.semispan)], 'k');

		plot(y_space, x_leading_edge(y_space), 'k')
		plot(y_space, x_trailing_edge(y_space), 'k')

		% Plot front and rear spar
		plot(y_space, x_front_spar(y_space), 'k')
		plot(y_space, x_rear_spar(y_space), 'k')

		grid on;
		axis equal;
	end

	total_volume = 0;
	num_stringers = starting_no_of_stringers; % This var will keep track of the number of stringers in the current panel;

	output.F_array = [];
	output.rib_array = [];
	while true
        
		bending_moment = bending_moment_dist(spanwise_station);

		box_width = geometry.box_width_func(spanwise_station);
		web_height = geometry.web_height_func(spanwise_station);

		stringer.flange_width = design_params.stringer_web_height*design_params.flange_to_web_ratio;
		stringer.cross_sec_area = design_params.stringer_thickness*design_params.stringer_web_height + 2*stringer.flange_width*design_params.stringer_thickness;% TODO: As appropriate for the selected stringer type

		K = 4 * (pi^2)/(12*(1-(material.v^2))); % From buckling of SS plate plot, taken for a/b -> infinity
		comp_load_per_length = bending_moment/(box_width*web_height);
		panel_thickness = ((comp_load_per_length*design_params.stringer_pitch^2)/(K*material.E))^(1/3);
		sigma_0 = comp_load_per_length/panel_thickness; % Critical Buckling Stress -> need to adjust this
		% for stringers using catchpole diagram

		[Kna,~,stringer_panel_area,eff_thickness] = stringer_panel_geometry(design_params.stringer_pitch,...
																			panel_thickness,...
																			design_params.stringer_thickness,...
																			design_params.stringer_thickness,...
																			stringer.flange_width,...
																			design_params.stringer_web_height,...
																			"Z_dh03",...
																			0);

		% eff_thickness = panel_thickness + stringer.cross_sec_area/design_params.stringer_pitch; 
		number_of_panels = num_stringers + 1;
		panel_eff_area = eff_thickness*design_params.stringer_pitch*number_of_panels;

		%% Catchpole diagram calculations
		t_s_over_t = design_params.stringer_thickness/panel_thickness;
		h_over_b = design_params.stringer_web_height/design_params.stringer_pitch;
		K_catchpole = catchpole_calculator(h_over_b, t_s_over_t, material.v,1);
		sigma_cr = K_catchpole/K*sigma_0;
		% sig_cr = K_catchpole*material.E*((panel_thickness/design_params.stringer_pitch)^2)
		% sig_0 = 3.62*material.E*((panel_thickness/design_params.stringer_pitch)^2)
		% sig_cr = (K_catchpole/K)*sig_0

		%equate sigma_cr to euler buckling to find optimum length

		rib_spacing = Kna*pi*sqrt(material.E/sigma_cr);
		F = sigma_cr * sqrt(rib_spacing/(comp_load_per_length*material.E));
		

		assert(F < 1, 'F stands for fuuucked')


		%% FARRAR efficiency factor
		A_s_over_bt = stringer.cross_sec_area/(design_params.stringer_pitch*panel_thickness);
		Fcalc = farrar_calculator(A_s_over_bt, t_s_over_t);
		rbs = comp_load_per_length*material.E*(Fcalc/sigma_cr)^2;

		% sg = Fcalc*sqrt((comp_load_per_length*material.E)/)
		
		%% Draw a rib
		spanwise_station = spanwise_station + rib_spacing;
                
		previous_stringers_to_cut = intercepts < spanwise_station & intercepts > spanwise_station - rib_spacing;

		if any(previous_stringers_to_cut) 
			intercepts(previous_stringers_to_cut) = (spanwise_station - rib_spacing); % cut the previous stringers which didnt reach far enough
			% now we have to redo the current iteration because the number of stringers changed (as some previous ones were cut)
			spanwise_station = spanwise_station - rib_spacing;
			num_stringers = num_stringers - sum(previous_stringers_to_cut, 'all');
            continue;
        end

        if spanwise_station > geometry.semispan
			total_volume = total_volume + panel_eff_area*(geometry.semispan - (spanwise_station - rib_spacing)); % panel untill the end
            break;
        end

		if doPlot
			plot([spanwise_station, spanwise_station], [x_leading_edge(spanwise_station), x_trailing_edge(spanwise_station)], 'r');
		end


		output.F_array = [output.F_array, F];
		output.rib_array = [output.rib_array, spanwise_station];
		
		% intercepts(stringers_to_cut) = spanwise_station;
		total_volume = total_volume + panel_eff_area*rib_spacing;
		total_volume = total_volume + geometry.A0*geometry.c(spanwise_station)^2*2E-3;	% account for ribs as well, might be kinda dodgy tho
																						% bcs we dont know rib thickness yet
	end

	if doPlot
		stringer_x_data = repmat(y_space, starting_no_of_stringers, 1)';
		stringer_y_data = stringer_func(y_space, stringer_x_space')';

		valid_indices = stringer_x_data < intercepts';
        stringer_x_data(~valid_indices) = NaN;
        stringer_y_data(~valid_indices) = NaN;
		plot(stringer_x_data, stringer_y_data, 'g');
		hold off;
	end

	output.total_volume = total_volume;
	output.total_weight = total_volume*material.rho;
end