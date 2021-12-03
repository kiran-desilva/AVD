clear

load('wing')

placer.nose_wheel.x = 0.5;

nose = linspace(0.3, 1.5, 10)
% nose = 0.3;
for n = nose
	placer.nose_wheel.x = n;
	% main_x = linspace(max(1, placer.nose_wheel.x), 11.736, 100)
	main_x = linspace(placer.nose_wheel.x, 11.736, 300)
	main_y = linspace(0, wing.b/2, 20)

	points = [];
	for x = main_x
		for y = main_y
			placer.main_wheel.x = x;
			placer.main_wheel.y = y;

			fail = 0;
			undercarriage;

			points = [points; x, y, fail];
		end
	end

	figure;
	scatter3(points(:, 2), points(:, 1), points(:, 3), [], points(:, 3));
	title(['Nose gear = ', num2str(n), ' m']);
	axis equal;
	caxis([0, 5]);
	points(:, 3)

end

function c = transform_arr(arr)
	c = [];
	for a = arr
		if a == 0
			c = [c; ]
		end
	end
end 