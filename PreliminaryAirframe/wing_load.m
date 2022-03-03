%% This function calculates shear, bending moment and torsion distributions for a given load case
%% SI units used throughout
%	n						-	load factor of the load case
%	Vinf_eas				-	equivalent airspeed 
%	wing_fuel_weight_kg		-	fuel weight in the wing
%	cl_dist					-	array representing discretised sectional lift coefficient distribution
%	spanwise_disc			-	array with spanwise points at which cl is sampled
%								THIS HAS TO BE UNIFORMLY SPACED!!!
function distributions = wing_load(n, Vinf_eas, wing_fuel_weight_kg, cl_dist, spanwise_disc, cm0, isLanding, doPlot, loadingCaseStr)
	load('loadcase.mat')

	% n = 1;
	% Vinf_eas = 100;
	% wing_fuel_weight_kg = 100;
	% cl_dist = ones(1, 1000);
	% spanwise_disc = linspace(0, 8.94/2, 1000);
	% isLanding = false;

	delta_s = spanwise_disc(2) - spanwise_disc(1);
	g = 9.81;
	rho0 = 1.225;
	Sref = 10.32;
	lbs_to_kg = 0.453592;

	spanwise_disc = spanwise_disc(1:end) - delta_s/2;
	assert(all(size(spanwise_disc) == size(cl_dist)), 'Got incompatible sizes for spanwise_disc and cl_dist');
	assert(all(diff(diff(spanwise_disc))) == 0, 'Spanwise_disc not a uniform distribution');
	assert(spanwise_disc(1) ~= 0, 'distribution wrong');

	syms x;
	c_root = 1.767;
	c_tip = 0.533;
	span = 8.94;
	semispan = span/2;
	t_c_ratio = 0.12;

	c(x) = -(c_root - c_tip)/semispan*x + c_root;
	% {
	%% Plot fuselage
	fus.r = 1.68/2;
	th = 0:pi/50:2*pi;
	xunit = fus.r * cos(th);
	yunit = fus.r * sin(th);
	color = 'k';

	plot(xunit, yunit, color);
	axis equal;
	grid on;


	%% Underbelly
	ubelly.x_at_r = 0.51;
	ubelly.y_at_r = -sqrt(fus.r^2 - ubelly.x_at_r^2);

	ubelly.x_at_c = 0;
	ubelly.y_at_c = -0.14-fus.r;

	syms A B
	f(x) = A*exp(15*x) + B;

	[solA, solB] = solve(f(ubelly.x_at_r) == ubelly.y_at_r, f(ubelly.x_at_c) == ubelly.y_at_c);
	f(x) = subs(f, [A B], [solA, solB]);
	x_range = linspace(ubelly.x_at_c, ubelly.x_at_r, 100);
	hold on;
	plot(x_range, f(x_range), color);
	plot(-x_range, f(x_range), color);
	hold off

	h(x) = t_c_ratio*c(x);

	x_heights = linspace(0, semispan, 100);
	top_edge(x) = h(x)/2 - fus.r;
	bottom_edge(x) = -h(x)/2 - fus.r;

	top_edge_intercept = double(vpasolve(top_edge == f, x, [0 semispan]));
	bottom_edge_intercept = double(vpasolve(bottom_edge == f, x, [0 semispan]));

	top_linspace = linspace(top_edge_intercept, semispan, 100);
	bottom_linspace = linspace(bottom_edge_intercept, semispan, 100);

	% hold on;
	% plot(top_linspace, top_edge(top_linspace), color);
	% plot(bottom_linspace, bottom_edge(bottom_linspace), color);
	% plot([semispan semispan], [bottom_edge(semispan) top_edge(semispan)], color)
	% hold off

	%% Wing aero loading
	% lift.gamma0 = 4*semispan*Vinf_eas*1.3/(pi*7.8);
	% lift.gamma(x) = lift.gamma0*sqrt(1 - (x/semispan).^2)
	% lift.L(x) = rho0*Vinf_eas*lift.gamma(x)
	% lift.L_dist = double(lift.L(spanwise_disc))
	lift.L_dist = 0.5*rho0*n*Vinf_eas^2*double(c(spanwise_disc)).*cl_dist;
	%double(trapz(spanwise_disc, lift.L_dist))

	%% Wing weight loading
	% M(x1, x2) = M_wing*V(x1, x2)/V_wing
	% V(x1, x2) = int_x1^x2{A0*c(x)^2*dx} where A0 is the non-dimensional area of the airfoil section
	syms x1 x2
	wing.A0 = 0.0655;
	wing.M = 189.6/2*lbs_to_kg;
	V(x1, x2) = int(wing.A0*c(x)^2, x, [x1, x2]);
	wing.volume = V(bottom_edge_intercept, semispan);
	wing.avg_density = wing.M/wing.volume;
	M(x1, x2) = wing.avg_density*V(x1, x2);
	x1_arr = spanwise_disc - delta_s/2;
	x2_arr = spanwise_disc + delta_s/2;
	wing.inertial_loading = -g*n*double(M(x1_arr, x2_arr))/delta_s;

	%% Undercarriage loading
	uc.spanwise_start = 1.18 - 0.049;
	uc.spanwise_end = 1.18 + 0.049 + 0.2;
	uc.wing_loading = 80*lbs_to_kg*g*n/abs(uc.spanwise_start - uc.spanwise_end);
	uc.loading(x) = piecewise((x >= uc.spanwise_start) & (x <= uc.spanwise_end), -uc.wing_loading, 0);

	%% Fuel loading
	% fuel.volume_in_wings_m3 = 0.5341;
	% fuel.fuel_weight = g*fuel.density*fuel.volume_in_wings_m3;
	fuel.density = 775; % kg/m^3
	fuel.full_wing_mass = fuel.density*V(bottom_edge_intercept, 0.95*semispan);
	fuel.percent_full = double(wing_fuel_weight_kg/fuel.full_wing_mass);
	fuel.loading(x1, x2) = piecewise((x1 >= bottom_edge_intercept) & (x2 <= 0.95*semispan), -g*n*fuel.density*fuel.percent_full*V(x1, x2)/delta_s, 0);

	%% Torsion dist
	% Positive torque clockwise
	% Needs to include:
	%	1.	Lift force at x_ac
	%	2.	Pitching moment coeff
	%	3.	Fuel force (assume at 0.5 x/c??) MIGHT BE BETTER TO ASSUME CoM of aerofoil
	%	4.	Wing weight at CoM of aerofoil
	%	5.	Undercarriage weight

	% Calculating shear force location
	x_front_spar_percent_c = 0.1;
	x_rear_spar_percent_c = 0.77;
	x_sc_assumption_percent_c = (0.1 + 0.77)/2;

	% Lift force calc
	x_ac_percent_c = 0.25;
	lift.torsional_load = -lift.L_dist.*(x_ac_percent_c - x_sc_assumption_percent_c).*c(spanwise_disc);

	% Pitching moment inclusion
	%c_mac = 1.2607;
	lift.pitching_moment_load = 0.5*rho0*Vinf_eas^2*c(spanwise_disc).^2*cm0;

	% Fuel tank contribution
	x_fuel_percent_c = 0.413575;
	fuel.torque(x1, x2) = -fuel.loading(x1, x2)*(x_fuel_percent_c - x_sc_assumption_percent_c)*(c(x1) + c(x2))/2;

	% Wing weight contribution
	wing.centroid_x_percent_c = 0.4136;
	wing.torsional_load = -wing.inertial_loading.*(wing.centroid_x_percent_c - x_sc_assumption_percent_c).*(c(x1_arr) + c(x2_arr))/2; % dimensionalise by avg cord between x1 and x2

	% UC weight contribution (update if isLanding)
	uc.attachment_point_percent_c = (c(1.18) - 0.1)/c(1.18);
	uc.torsional_load(x) = -uc.loading(x)*(uc.attachment_point_percent_c - x_sc_assumption_percent_c)*c(x);

	if (isLanding)
		xcg = 4.55; % measured from the nose in m
		xG = 5.16;
		xT = 10.3474;
		W = 3128.2*0.85; % landing weight is 85% of mass takeoff weight
		F_uc = W*g*(xcg - xT)/(2*(xG - xT)); % Points up!! => positive

		uc.loading(x) = piecewise((x >= uc.spanwise_start) & (x <= uc.spanwise_end), F_uc/abs(uc.spanwise_start - uc.spanwise_end), 0);
		uc.torsional_load(x) = -uc.loading(x)*(uc.attachment_point_percent_c - x_sc_assumption_percent_c)*c(x);

		lift.L_dist = 0;
		lift.torsional_load = 0;

		lift.pitching_moment_load = 0;
    end
    
    combined_load = double(uc.loading(spanwise_disc) + fuel.loading(x1_arr, x2_arr) + lift.L_dist + wing.inertial_loading);
    combined_force = movsum(combined_load, 2)*delta_s/2;
    combined_force = [combined_force(2:end), 0];
    
	shear_dist = sum(combined_force) - cumsum(combined_force) + combined_force;

	temp = movsum(shear_dist, 2)*delta_s/2;
	dM = [temp(2:end), 0];

	bm_dist = sum(dM) - cumsum(dM) + dM;
    
    %lift.pitching_moment_load = 0;
	torsional_dist = double(uc.torsional_load(spanwise_disc) + wing.torsional_load + lift.pitching_moment_load + lift.torsional_load + fuel.torque(x1_arr, x2_arr));
    
    temp = movsum(torsional_dist, 2)*delta_s/2;
    dT = [temp(2:end), 0];
	torque_dist = sum(dT) - cumsum(dT) + dT;

	distributions.points = spanwise_disc;
	distributions.bm = bm_dist;
	distributions.torque = torque_dist;
	distributions.shear = shear_dist;

	if doPlot
		figure;
		subplot(5,1,1)
		hold on;
		plot(spanwise_disc, torsional_dist);
		grid on;
		title(append(loadingCaseStr, "Combined Torque"))
		hold off;

		subplot(5,1,2)
		hold on;
		plot(spanwise_disc, torque_dist);
		grid on;
		title(append(loadingCaseStr, "Torque Dist"))
		hold off;

		subplot(5,1,3)
		hold on;
		plot(spanwise_disc, combined_force);
		% xline(bottom_edge_intercept, 'r')
		% plot(X_hatch, Y_hatch, 'r')
		xlim([0 semispan]);
		% ylim([0 max(combined_loading)*1.1])
		grid on;
		title(append(loadingCaseStr, "Combined Loading"))
		hold off

		subplot(5,1,4)
		hold on;
		plot(spanwise_disc, shear_dist);
		% xline(bottom_edge_intercept, 'r')
		grid on;
		title(append(loadingCaseStr, "Shear Loading"))
		xlim([0 semispan]);
		% ylim([0 max(shear_dist)*1.1])
		hold off;

		subplot(5,1,5);
		hold on;
		plot(spanwise_disc, bm_dist);
		% xline(bottom_edge_intercept, 'r')
		grid on;
		title(append(loadingCaseStr, "Bending Moment"))
		xlim([0 semispan]);
		% ylim([0 max(bm_dist)*1.1])
		hold off;

		figure;
		subplot(2,1,1);
		hold on;
		plot(spanwise_disc, lift.L_dist);
		% xline(bottom_edge_intercept, 'r')
		grid on;
		title(append(loadingCaseStr, "Lift Distribution"))
		xlim([0 semispan]);
		% if (~isLanding)
		%     ylim([0 max(lift.L_dist)*1.1]);
		% end
		hold off;

		subplot(2,1,2);
		hold on;
		plot(spanwise_disc, wing.inertial_loading);
		% xline(bottom_edge_intercept, 'r')
		grid on;
		title(append(loadingCaseStr, "Wing Weight Distribution (Debug only)"))
		xlim([0 semispan]);
		% ylim([max(double(wing.inertial_loading))*1.1, 0])
		hold off;

		% hold on;
		% plot(bottom_linspace, fuel.loading(bottom_linspace) + wing.loading(bottom_linspace) + uc.loading(bottom_linspace) - fus.r);
		% hold off;
	end

end