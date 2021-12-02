clear

load('wing')

placer.nose_wheel.x = 0.5;

nose = linspace(0, 2, 10)
main_x = linspace(max(1, placer.nose_wheel.x), 11.736, 500)
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

scatter3(points(:, 2), points(:, 1), points(:, 3), [], points(:, 3));
axis equal;
caxis([0, 5]);
points(:, 3)

function c = transform_arr(arr)
	c = [];
	for a = arr
		if a == 0
			c = [c; ]
		end
	end
end 