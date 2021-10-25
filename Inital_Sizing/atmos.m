function ret = atmos(height, idx) 
	[a, b, c, d] = atmosisa(height);
	mat = [a, b, c, d];
	ret = mat(idx);
end