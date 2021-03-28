function [Rbar, Phi, F, Phi_Phi, Phi_F,Phi_R, A_e, B_e, C_e] = mpcPojacanjaMIMO(Ap,Bp,Cp,Nc,Np,r_omega)
%mpcPojacanjaMIMO Pronalazi pojacanja  MPC regulatora za zadati objekat p-plant
%
% [Phi_Phi, Phi_F,Phi_R, A_e, B_e, C_e] = mpcgain(Ap,Bp,Cp,Nc,Np)
%
% Ap, Bp, Cp - Diskretne matrice stanja objekta
% A_e, B_e, C_e - Izmenjene, poboljsanje (Augmented) matrice 
% Nc - Horizont upravljanja
% Np - Horizont predikcije
% r_omega - Tezinski koeficient optimizacije upravljackog signala
% Optimalno dU => dU = inv(Phi_Phi)*(Phi_R-Phi_F*x(ki))
% gde je x(ki) - vektor stanja u trenutnu ki

%*****************************************************************
% Autor: M.Lomovic
%*****************************************************************
[m1,n1] = size(Cp);
[n1,n_in] = size(Bp);

%Kreiranej izmenjenih matrica modela u odnosu na odstupanja
A_e = eye(n1+m1,n1+m1);
A_e(1:n1,1:n1) = Ap;
A_e(n1+1:n1+m1,1:n1) = Cp*Ap;

B_e = zeros(n1+m1,n_in);
B_e(1:n1,:) = Bp;
B_e(n1+1:n1+m1,:) = Cp*Bp;

C_e = zeros(m1,n1+m1);
C_e(:,n1+1:n1+m1) = eye(m1,m1);

[m1,n1] = size(C_e);
[n1,n_in] = size(B_e);

n = n1+m1;

% Kreiranje matrice F = [CA;CA^2;...;CA^Np]
F = zeros(m1*Np,n1);
stepen = 1;
for i = 1:m1:Np*m1   
    F(i:i+m1-1,:) = C_e*A_e^stepen;
    stepen = stepen +1;
end

% Kreiranje matrice Phi koja je Toeplitz-ova matrica
h = zeros(m1*Nc,n1);
%h(1:m1,:) = C_e
stepen = 0;
for i = 1:m1:Np*m1   
    h(i:i+m1-1,:) = C_e*A_e^stepen;
    stepen = stepen +1;
end

v = h*B_e;
Phi = zeros(Np*m1,Nc*m1); % dimezija Phi
Phi(:,1:m1) = v(:,1:m1); %Prva kolona Phi


for i = m1+1:n_in:Nc*n_in
    Phi(:,i:i+m1-1)=[zeros(i-1,m1);v(1:Np*m1-i+1,:)]; %Toeplitz-ova matrica
end



% Kreiranje ostalih matrica za predikciju
%BarRs - Rs nadvuceno
BarRs = zeros(Np*m1,n_in);
for i = 1:m1:Np*m1
    BarRs(i:i+m1-1,:) = eye(n_in);
end
Phi_Phi = Phi'*Phi;
Phi_F = Phi'*F;  
Phi_R = Phi'*BarRs;
%Rbar = zeros(Np*m1,Nc*m1);
Rbar = r_omega*eye(Nc*m1,Nc*m1);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Provera
%Fmoje = [C_e*A_e; C_e*A_e^2; C_e*A_e^3; C_e*A_e^4; C_e*A_e^5;]

%o = [0 0;0 0];

%Fi = [   C_e*B_e          o                  o                  o               o;
%      C_e*A_e^1*B_e     C_e*B_e             o                  o                o;
%      C_e*A_e^2*B_e     C_e*A_e^1*B_e      C_e*B_e             o                o;
%      C_e*A_e^3*B_e     C_e*A_e^2*B_e      C_e*A_e^1*B_e     C_e*B_e            o;
%      C_e*A_e^4*B_e     C_e*A_e^3*B_e      C_e*A_e^2*B_e     C_e*A_e^1*B_e     C_e*B_e];














