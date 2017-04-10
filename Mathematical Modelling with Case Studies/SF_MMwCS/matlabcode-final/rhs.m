function udot = rhs(t, u)
global a1 a2;

R=u(1); B=u(2); 
Rdash = -a1*B;
Bdash = -a2*R;
udot = [Rdash; Bdash];

