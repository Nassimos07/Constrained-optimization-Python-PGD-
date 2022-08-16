%Partie_2 :cas où les paramètres a, b, c et d ne sont pas constants 
%Initialisation 
a=1;
b=0.7;
c=2;
d=-1;

N=100 ;
U=randi([0 1],1,N);  
% Y=zeros(1,N);
Y(1)=U(1);
Y(2)=U(2);


for i=3:N;
    a_k= a*(1-0.5*exp(-i));
    b_k= b*(1-0.5*exp(-i));
    c_k= c*(1-0.5*exp(-i));
    d_k= d*(1-0.5*exp(-i));

    Y(i)=c_k*U(i-1)+d_k*U(i-2)-a_k*Y(i-1)-b_k*Y(i-2);% creation de la sortie 

end; %

bruit_faible_amplitude= 0.1*randi([0 1],1,N); 

Y=Y+bruit_faible_amplitude; %Ajout de la bruit à la sortie
% Identification par moindres carrés recursifs
%initialisation 
alpha=1000;
% initialiser teta par les premier valeur 
A=[-Y(5) -Y(4) U(5) U(4);-Y(4) -Y(3) U(4) U(3);-Y(3) -Y(2) U(3) U(2);-Y(2) -Y(1) U(2) U(1)];
B=[Y(6);Y(5);Y(4);Y(3)];
teta=pinv(A)*B;
P=alpha*eye(4);
TETA=teta;
% pour le cas ou a,b,c et d sont variables on utilise l'algorithme avec
% facteur d'oublie
%odre d'itiration = nombres des échantillons 


lambda = .95;  % facteur d'oublie 

for k=7:N;  
     h=[-Y(k-1) -Y(k-2) U(k-1) U(k-1)]';
     G    =         P*h*inv(lambda+h'*P*h);
     teta =         teta+G*(Y(k)-h'*teta);
     P    =         (1/lambda)*(eye(4)-G*h')*P;
     TETA =         [TETA,teta];
end;
a_estimated=teta(1)
b_estimated=teta(2)
c_estimated=teta(3)
d_estimated=teta(4)
display("\n")
Erreur=[a;b;c;d]-teta

subplot(221)
plot(1:N-5,TETA(1,:),'r')
subplot(222)
plot(1:N-5, TETA(2,:),'b')
subplot(223)
plot(1:N-5,TETA(3,:),'g')
subplot(224)
plot(1:N-5, TETA(4,:),'y')