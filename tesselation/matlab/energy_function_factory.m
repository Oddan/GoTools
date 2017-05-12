function efun = energy_function_factory(type, radius)

   switch(type)
     case 'simple'
       efun = @(l) simple_energy_function(l, radius);
     otherwise 
       error('unknown energy function type.');
   end
   
end

function [e, de, de2] = simple_energy_function(l, radius)
   
   tmp = max(radius - l, 0);
   
   e = tmp.*tmp;
   de = -2 * tmp;
   de2 = 2 * (tmp > 0);
   
end
