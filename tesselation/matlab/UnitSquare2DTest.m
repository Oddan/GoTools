function hist = UnitSquare2DTest(radius)
   
  P = [0.6, 0.55; 0.4 0.4; 0.9, 0.1];
  corners = [ 0 0; 1 0; 1 1; 0 1];
  
  hist = struct('P', P, 'err', []);
  small = 0.001; %0.5;
  MAX_ITER = 50;
  threshold = 1e-12;
  
  for i = 1:MAX_ITER
     
     [e, d, dd] = energy(radius, P, corners);

     error = norm(d, inf);
     hist.err = [hist.err; error];

     if i < 10
        update = -d * small;
     else
        % eliminate points that do not contribute to energy
        ix = ((abs(d)<threshold) & (sum(abs(dd), 2) < threshold));
        
        if sum(ix)  == 0 
           update = -dd\d;
        else
           keep = ~ix;
           tmp = -dd(keep, keep)\d(keep);
           update = 0 * d;
           update(keep) = tmp;
        end
     end
     update = reshape(update, [], 2);

     
     
     P_old = P;
     P = P + update;

     % ix_outside = (P<= 0) | (P >= 1);
     % P(ix_outside) = P_old(ix_outside);
     
     backtrack = update;
     count = 0;
     while any((P(:) <= 0) | (P(:) >= 1))
        backtrack = backtrack / 2;
        P = P - backtrack;
        count = count + 1;
        if count == 1000
           keyboard;
        end
     end
     
     hist.P = [hist.P; P];
     
     fprintf('Iteration %i - error: %d \n', i, error);
     
     if error < threshold
        break;
     end
  end

  P
         
end



