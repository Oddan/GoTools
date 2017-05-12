%polygon = [0,0; 1, 0; 1, 1; 0, 1];
polygon = [0.5,0; 1, 0.5; 0.5, 1; 0, 0.5];

points = [0.5, 0.5];

efun = energy_function_factory('simple', 1);
dfun = @euclidian_distance;

[E, dE] = polygon_energy(polygon, points, dfun, efun);

%mirror_point_energy(seg, points, @euclidian_distance, efun);

