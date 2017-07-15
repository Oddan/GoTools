function plot_tri(pts, tri)

   clf; hold on
   
   for i = 1:size(tri, 1)
      
      ixs = [tri(i,:), tri(i,1)]';
      
      plot(pts(ixs,1), pts(ixs,2), '*-r');
      
   end
   
end

