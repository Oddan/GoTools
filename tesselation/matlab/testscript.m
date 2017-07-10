polygon = [0,0; 1, 0; 1, 1; 0, 1];
%polygon = [0.5,0; 1, 0.5; 0.5, 1; 0, 0.5];
polygon = [0.5, -0.5; 1.5, 0.5; 0.5, 1.5; -0.5, 0.5]; % rotated square centered on (0.5, 0.5)
polygon = [0.5, -0.5; 1.5, 0.5; 1.5, 1.5; 0.5, 1.5; -0.5, 0.5]; % rotated house


points = [0.5, 0.5];
points = [0.1, 0.1; 0.5, 0.5];
          

% Really, the radius should be chosen by dividing total polygon area by
% number of points, computing the equivalent radius (divide by pi and take
% square root), and then multiplying by a "reasonable small" value
% (perhaps 4 or 5) which represents the average number of points in
% neighborhood.
% A larger value of the radius will require more computational time, and
% force points closer to boundary.
R = sqrt(polygon_area(polygon) / (numel(points)/2)/pi) * 5;
efun = energy_function_factory('simple', R); 
%efun = energy_function_factory('steep', 1);
dfun = @euclidian_distance;

[E, dE] = polygon_energy(polygon, points, dfun, efun);

x = cgrad2D(@(x) polygon_energy(polygon, reshape(x(:), [], 2), dfun, efun), points(:), polygon);


%mirror_point_energy(seg, points, @euclidian_distance, efun);

