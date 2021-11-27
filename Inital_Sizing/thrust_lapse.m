function beta = thrust_lapse(height_m, mach_number, bypass_ratio)
	[~, ~, ~, rho0] = atmosisa(0);
	[T, a, P, rho] = atmosisa(height_m);


	% From Errikos' notes
	beta_altitude = 1.438*rho./rho0;

	% From https://www.fzt.haw-hamburg.de/pers/Scholz/arbeiten/TextSchulzDipl.pdf page 22 (Ultimately references Torenbeek)
	G = 0.9; % For low BPR engine

	beta_mach = 1 - (0.45*mach_number*(1 + bypass_ratio))./(sqrt(G*(1 + 0.75*bypass_ratio))) + mach_number.^2*(0.6 + 0.11*bypass_ratio/G);

	beta = beta_altitude .* beta_mach;

end