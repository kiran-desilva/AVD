function beta = thrust_lapse(height_ft, mach_number, bypass_ratio)
	height_m = distdim(height_ft, 'ft', 'm');
	[~, ~, ~, rho0] = atmosisa(0);
	[T, a, P, rho] = atmosisa(height_m);


	% From Errikos' notes
	if height_m > 11000
		beta_altitude = 1.439*rho./rho0;
	else
		beta_altitude = (rho./rho0)^0.7;
	end

	% From https://www.fzt.haw-hamburg.de/pers/Scholz/arbeiten/TextSchulzDipl.pdf page 22 (Ultimately references Torenbeek)
	G = 0.9; % For low BPR engine

	beta_mach = 1 - (0.45*mach_number*(1 + bypass_ratio))./(sqrt(G*(1 + 0.75*bypass_ratio))) + mach_number.^2*(0.6 + 0.11*bypass_ratio/G);

	beta = beta_altitude .* beta_mach;

end