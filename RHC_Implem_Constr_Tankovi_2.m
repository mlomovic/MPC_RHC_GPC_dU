%Modeliranje objekta
%
% A1 - Površina dna posude 1 [m^2]
% R1 - Otpor isticanja preko ventila V12 izme?u dva tanka [s/m^2]
% h1 - trenutna visina vodenog stuba u posudi 1 [m]
% A2 - Površina dna posude 1 [m^2]
% R2 - Otpor isticanja preko ventila Vi posle taanka 2 [s/m^2]
% h2 - trenutna visina vodenog stuba u posudi 2 [m]
%
%  dh1       -1             1          1
%  --- =   -------*h1 +  -------*h2 + ---*Qu1
%  dt       A1*R1          A1*R1       A1
%
%  dh2        1             1       1      1        1
%  --- =   -------*h1 - ( ----- + ----- )----*h2 + ---*Qu2
%  dt       A2*R1           R1      R2    A2        A2
%
% Model objekta
%                                       
% dx/dt = Ax(t) + Bu(t)
%  y(t) = Cx(t) + Du(t)
%
clear, clc;

%Fizicke velicine
A1 = 1;           %Povrsina dna tanka 1 [m^2]
R1 = 1;          %Otpor isticanja preko ventila 12 [s/m^2]
A2 = 1;           %Povrsina dna tanka 2 [m^2]
R2 = 1;          %Otpor isticanja preko ventila 12 [s/m^2]

dt = 0.1;       %Perioda odabiranja diskretnog sistema [s]



%Kontinualni sistem
 
Ac = [-1/(A1*R1) 1/(A1*R1);
     1/(A2*R1)  -(1/A2*R1 + 1/A2*R2)];
 
Bc = [1/A1 0;
     0 1/A2];
 
Cc = [1 0;
     0 1];

Dc = [0 0;
     0 0];
 
 
sis = ss(Ac,Bc,Cc,Dc)

sisd = c2d(sis,dt,'tustin')
[Ap,Bp,Cp,Dp] = ssdata(sisd)
 
Np=10;            %Horizont predikcije stanja (izlaza)
Nc=10;            %Horizont predikcije upravljanja Nc <= Np
r_omega=0.01;    %Tezinski koeficient upravljanja (r_omega >= 0 za svako ki)
N_sim=15;        %Duzina simulacije (broj trenutaka odabiranja) N_sim >= Np

% Ogranicenja
u_max = [4; 4];
u_min = [-3;-3];
deltau_max = [-1.5; -1.5];
deltau_min = [1.5; 1.5];


[R_bar,Phi,F,Phi_Phi,Phi_F,Phi_R,A_e,B_e,C_e] = mpcPojacanjaMIMO(Ap,Bp,Cp,Nc,Np,r_omega);

Kmpc = inv(Phi_Phi+R_bar)*Phi_F; % Pojacanje MPC-a
Ky = inv(Phi_Phi+R_bar)*Phi_R;


[m,n] = size(Cp);
[n,n_in]=size(B_e);

% Zadavanje zeljenih vrednosti izlaza
r1=ones(N_sim,1);
%r1(50:150,1)=1;
r2=ones(N_sim,1)*0.4;
%r2(200:300,1)=1;
r=[r1 r2];


% Pocetni uslovi
xm=[0;0];
Xf=zeros(n,1);

u=[0;0]; % u(k-1) =0
y=zeros(m,1);

% Matrice za  optimizaciju sa ogranicenjima
A_cons = zeros(4*Nc*n_in,Nc*n_in);
for i=1:Nc*n_in
    for j=1:Nc*n_in 
        if i==j
            A_cons(i,j)=1;   
            A_cons(i+Nc*n_in,j)=-1;
        end;
        if i>=j 
            A_cons(i+2*Nc*n_in,j)=1;
            A_cons(i+3*Nc*n_in,j)=-1;
        end;        
    end;
end;

b = zeros(4*Nc*n_in,1);

for i = 1 : n_in : Nc*n_in
    b(i:i+n_in-1,1) = deltau_min;
end
for i = Nc*n_in+1 : n_in : 2*Nc*n_in-1
     b(i:i+n_in-1,1) = -deltau_max;
end
for i = 2*Nc*n_in+1 : n_in : 3*Nc*n_in-1
     b(i:i+n_in-1,1) = u_max-u;
end
for i = 3*Nc*n_in+1 : n_in : 4*Nc*n_in-1
     b(i:i+n_in-1,1) = -u_min+u;
end

%b1=[deltau_min; deltau_min; deltau_min; -deltau_max; -deltau_max; -deltau_max; u_max-u; u_max-u; u_max-u; -u_min+u; -u_min+u; -u_min+u];

for kk=1:N_sim;
DeltaU1 = inv(Phi_Phi+R_bar) * (Phi_R*r(kk,:)'-Phi_F*Xf);
H = Phi_Phi+R_bar;
f = -Phi_R*r(kk,:)'+Phi_F*Xf;

%b1=[deltau_min; deltau_min; deltau_min; -deltau_max; -deltau_max; -deltau_max; u_max-u; u_max-u; u_max-u; -u_min+u; -u_min+u; -u_min+u];

for i = 2*Nc*n_in+1 : n_in : 3*Nc*n_in-1
     b(i:i+n_in-1,1) = u_max-u;
end
for i = 3*Nc*n_in+1 : n_in : 4*Nc*n_in-1
     b(i:i+n_in-1,1) = -u_min+u;
end

DeltaU = QPhild(H,f,A_cons,b);
deltau=[eye(n_in),zeros(n_in,Nc*n_in-n_in)]*DeltaU;%(1,1);
u=u+deltau;
u1(kk,:)=u;
y1(kk,:)=y;
xm_old=xm;
xm=Ap*xm+Bp*u;
y=Cp*xm;
Xf=[xm-xm_old;y];
deltau_crtanje(kk,1:2)=deltau';
end


k=0:(N_sim-1);
figure
s(1) = subplot(311)
plot(k,y1(:,1),'-b',k,y1(:,2),'-r',k,r1,'--b',k,r2,'--r'),grid;
title(s(1),['Upravljanje MPC-om sa parametrima: N_c = ',num2str(Nc),', N_p = ',num2str(Np),', r_\omega = ',num2str(r_omega), ', T = ',num2str(dt),'[s]']);
xlabel('Trenutak odabiranja [k]')
ylabel('Izlaz')
legend('y_1','y_2','r_1','r_2')
s(2) = subplot(312)
stairs(k,u1),grid;
xlabel('Trenutak odabiranja [k]')
ylabel('Upravljanje')
legend('u_1','u_2')
s(3) = subplot(313)
plot(k,deltau_crtanje),grid;
xlabel('Trenutak odabiranja [k]')
ylabel('Prirastaj upravljanja')
legend('\Deltau_1', '\Deltau_2')



[K,S,E] = dlqr(Ap,Bp,Cp,[0.1 0; 0 0.1])

sisdc = ss((Ap-Bp*K),Bp,Cp,Dp,0.1)
%figure
%step(sisdc),grid

