function [Phi, F, Phi_Phi, Phi_F,Phi_R, A_e, B_e, C_e] = mpcPojacanja(Ap,Bp,Cp,Nc,Np)
%MPCGAIN Pronalazi pojacanja  MPC regulatora za zadati objekat p-plant
%
% [Phi_Phi, Phi_F,Phi_R, A_e, B_e, C_e] = mpcgain(Ap,Bp,Cp,Nc,Np)
%
% Ap, Bp, Cp - Diskretne matrice stanja objekta
% A_e, B_e, C_e - Izmenjene, poboljsanje (Augmented) matrice 
% Nc - Horizont upravljanja
% Np - Horizont predikcije
%
% Optimalno dU => dU = inv(Phi_Phi)*(Phi_R-Phi_F*x(ki))
% gde je x(ki) - vektor stanja u trenutnu ki

%*****************************************************************
% Autor: M.Lomovic
%*****************************************************************


[m1,n1] = size(Cp);
[n1,n_in] = size(Bp);

A_e = eye(n1+m1,n1+m1);
A_e(1:n1,1:n1) = Ap;
A_e(n1+1:n1+m1,1:n1) = Cp*Ap;

B_e = zeros(n1+m1,n_in);
B_e(1:n1,:) = Bp;
B_e(n1+1:n1+m1,:) = Cp*Bp;

C_e = zeros(m1,n1+m1);
C_e(:,n1+1:n1+m1) = eye(m1,m1);

n = n1+m1;
h(1,:) = C_e;
F(1,:) = C_e * A_e;

for kk = 2:Np
    h(kk,:) = h(kk-1,:)*A_e;
    F(kk,:) = F(kk-1,:)*A_e;
end

v = h*B_e;
Phi = zeros(Np,Nc); % dimezija Phi
Phi(:,1) = v; %Prva kolona Phi

for i = 2:Nc
    Phi(:,i)=[zeros(i-1,1);v(1:Np-i+1,1)]; %Toeplitz-ova matrica
end

BarRs = ones(Np,1);
Phi_Phi = Phi'*Phi;
Phi_F = Phi'*F;
Phi_R = Phi'*BarRs;
