function [val, der, der2] = energy(radius, points, boundary_points) 

   N = size(points, 1); % number of interior points
   all_points = [points; boundary_points]; 
   
   [X, Y] = ndgrid(all_points(:,1), all_points(:,2));
   Y = Y'; 
   
   dX = X - X'; % x-differences between points
   dY = Y - Y'; % y-differences between points

   D = sqrt(dX.^2 + dY.^2); % entry (i,j) : distance between point p_i and p_j
   
   [evals, devals, ddevals] = dist_energy_2D(D, radius);
   
   val = 0.5 * sum(evals(:));
   
   Dm = D + eye(size(D,1)); % modified D to avoid 0/0 in expressions below
   Edx = devals .* dX ./ Dm; 
   Edy = devals .* dY ./ Dm;
   
   der = [sum(Edx, 2), sum(Edy, 2)];
   der = der(1:N, :);
   der = der(:);

   dxx = -1./(Dm.^2) .* (ddevals .* (dX.^2) + devals .* (D - (dX.^2)./Dm));
   dxx(logical(eye(size(dxx,1)))) = -sum(dxx, 2);
   dxx = dxx(1:N, 1:N);
   
   dyy = -1./(Dm.^2) .* (ddevals .* (dY.^2) + devals .* (D - (dY.^2)./Dm));
   dyy(logical(eye(size(dyy, 1)))) = -sum(dyy, 2);
   dyy = dyy(1:N, 1:N);
   
   dxy = dX .* dY ./ (Dm.^2) .* (devals./Dm - ddevals);
   dxy(logical(eye(size(dxy, 1)))) = -sum(dxy, 2);
   dxy = dxy(1:N, 1:N);

   der2 = [dxx, dxy; dxy, dyy];
           
end

function [val, der, der2] = dist_energy_2D(dist, radius)

   r = radius - dist;
   r(logical(eye(size(r)))) = 0; % no energy associated with a point and itself
   
   val  = r.*r;
   der  = -2 * r;
   der2 = 0*r+2;
   
   val (r<=0) = 0;
   der (r<=0) = 0;
   der2(r<=0) = 0;
   
end

