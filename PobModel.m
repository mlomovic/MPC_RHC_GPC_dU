% Kontinualni model u prostoru stanja
Ac = [0 1 0;
      3 0 1;
      0 1 0];
  
Bc = [1 1 3]' ;

Cc = [0 1 0];

Dc = zeros(1,1);


%Diskretizovani model u prostoru stanja c2dm - 'zoh' 

dt = 1; %perioda odabiranja
[Ad,Bd,Cd,Dd] = c2dm(Ac,Bc,Cc,Dc,dt);
[A, B, C] = FunPobModel(Ad, Bd, Cd)

%Odredjivanje dimenzija matrica sistema radi nalazenja broja ulaza i izlaza
%u sistem

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
