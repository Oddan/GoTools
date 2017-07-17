function plot_tri(pts, tri)

   clf; hold on
   
   for i = 1:size(tri, 1)
      
      ixs = [tri(i,:), tri(i,1)]';
      
      if size(pts,2) == 2
         plot(pts(ixs,1), pts(ixs,2), '*-r');
      else
         plot3(pts(ixs,1), pts(ixs,2), pts(ixs, 3), '*-r');
      end
   end
   
end

