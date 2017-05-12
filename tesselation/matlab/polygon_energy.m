function [E, dE] = polygon_energy(poly, points, dfun, efun)
% @@ assuming convex polygon, for now

   %% counting points and edges
   K = size(poly, 1);    % number of edges
      
   %% compute internal energy
   [E, dE] = interpoint_energy(points, points, dfun, efun, 'internal');
   
   %% adding boundary energy
   
   % boundary polygon energy
   [BE, dBE] = interpoint_energy(points, poly, dfun, efun, 'boundary');
   
   BOUNDARY_WEIGHT = 2; % @@ chosen because otherwise points tend to move too
                        % close to boundary points.  Perhaps this should be a
                        % tuneable parameter
   E  = E  + BOUNDARY_WEIGHT * BE; 
   dE = dE + BOUNDARY_WEIGHT * dBE;
   
   % mirror energy
   segments = [poly; poly(1,:)]; % closing polygon
   for i = 1:K
      cur_segment = segments(i:i+1,:);
  
      mpoints = mirror_points(cur_segment, points);
      [ME, dME] = interpoint_energy(points, mpoints, dfun, efun);
      
      E = E + ME;
      dE = dE + dME;
   end
end
