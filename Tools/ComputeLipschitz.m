function L = ComputeLipschitz(alpha,beta,eta,p,q,N)


L1 = p*alpha^(p-2) / beta^p;
L2 = p/(2*alpha^2) * max(1, (N*alpha^p/beta^p)^2);
L3 = (q-1)/eta^2;

L = L1 + L2 + L3;