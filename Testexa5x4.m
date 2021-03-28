clear, clc;

%First row
n11=0.34;
d11=[0.85 1];
n12=0.21;
d12=[0.42 1];
n13=0.50*[0.50 1];
d13=[12 0.4 1];
n14=0;
d14=1;
n15=6.46*[0.9 1];
d15=[0.07 0.3 1];
%second row
n21=-0.41;
d21=[2.41 1];
n22=0.66;
d22=[1.51 1];
n23=-0.3;
d23=[1.45 1];
n24=0;
d24=1;
n25=-3.72;
d25=[0.8 1];
%third row
n31=0.3;
d31=[2.54 1];
n32=0.49;
d32=[1.54 1];
n33=-0.71;
d33=[1.35 1];
n34=-0.20;
d34=[2.71 1];
n35=-4.71;
d35=[0.008 0.41 1];
%fourth row
n41=0;
d41=1;
n42=0;
d42=1;
n43=0;
d43=1;
n44=0;
d44=1;
n45=1.02;
d45=[0.07 0.31 1];



dt=1; %sampling interval
Gs=tf({n11 n12 n13 n14 n15;n21 n22 n23 n24 n25;...
n31 n32 n33 n34 n35;n41 n42 n43 n44 n45},...
{d11 d12 d13 d14 d15;d21 d22 d23 d24 d25;...
d31 d32 d33 d34 d35;d41 d42 d43 d44 d45});
Gsmin=ss(Gs,'min');
[Ac,Bc,Cc,Dc]=ssdata(Gsmin);
[Ap,Bp,Cp,Dp]=c2dm(Ac,Bc,Cc,Dc,dt,'zoh');
[m1,n1]=size(Cp);
[n1,n_in]=size(Bp);

a1=0.5;
a2=0.5;
a3=0.5;
a4=0.5;
a5=0.5;
N1=10;
N2=10;
N3=10;
N4=10;
N5=10;
a=[a1 a2 a3 a4 a5];
N=[N1 N2 N3 N4 N5];
Np=100;

%%%%%%%%%%%%%%%%
%Augment state equations
%%%%%%%%%%%%%%%%
A_e=eye(n1+m1,n1+m1);
A_e(1:n1,1:n1)=Ap;
A_e(n1+1:n1+m1,1:n1)=Cp*Ap;
B_e=zeros(n1+m1,n_in);
B_e(1:n1,:)=Bp;
B_e(n1+1:n1+m1,:)=Cp*Bp;
C_e=zeros(m1,n1+m1);
C_e(:,n1+1:n1+m1)=eye(m1,m1);
Q=C_e'*C_e;
R=0.1*eye(n_in,n_in);

[Omega,Psi]=dmpc(A_e,B_e,a,N,Np,Q,R);
L_m=zeros(n_in,sum(N));
[A1,L0]=lagd(a(1),N(1));
L_m(1,1:N(1))=L0';
In_s=1;
for jj=2:n_in;
[Al,L0]=lagd(a(jj),N(jj));
In_s=N(jj-1)+In_s;
In_e=In_s+N(jj)-1;
L_m(jj,In_s:In_e)=L0';
end
K=L_m*(Omega\Psi);
Acl=A_e-B_e*K;


[X,Y,Z]=dlqr(A_e,B_e,Q,R);
figure
plot(Z,'ro')
hold on
plot(eig(Acl),'b*')




























