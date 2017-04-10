de1 = 'beta*X - c1*X*Y - d1*X^2= 0';
de2 = 'c2*X*Y - alpha*Y -d2*Y^2= 0';
soln = solve(de1, de2, 'X', 'Y');
disp(soln.X); disp(soln.Y);
