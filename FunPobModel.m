function [A, B, C] = FunPobModel(Ad, Bd, Cd)

[m1,n1] = size(Cd);
[n1,n_ni] = size(Bd);


%Kreiranje poboljsanog (Augmented) modela A_e, B_e, C_e, D_e

A_e = eye(n1+m1,n1+m1);
A_e(1:n1,1:n1) = Ad;
A_e(n1+1:n1+m1,1:n1) = Cd*Ad;

B_e = zeros(n1+m1,n_ni);
B_e(1:n1,:) = Bd;
B_e(n1+1:n1+m1,:) = Cd*Bd;

C_e = zeros(m1,n1+m1);
C_e(:,n1+1:n1+m1) = eye(m1,m1);

A=A_e;
B=B_e;
C=C_e;

end

