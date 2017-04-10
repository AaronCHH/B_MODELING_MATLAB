function dydt = SP(t,y,k1,k2,DOsat,fBOD)
  k3 = k1*(1 + 0.5*sin(t*(pi+pi)));
  dydt(1) = fBOD-k3*y(1);
  dydt(2) = k2*(DOsat-y(2))-k3*y(1);
end