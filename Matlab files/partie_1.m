% On considère les paramètres a=1, b=0.7, c=2 et d=-1

%Initialisation 
display("\nles parmetres sont\n")
a=1
b=0.7
c=2
d=-1


N=100 ;
%-----------------------Question 1 -----------------
U=randi([0 1],1,N);  
%-----------------------Question 2 -----------------
% Y=zeros(1,N);
Y(1)=U(1);
Y(2)=U(2);

for i=3:N;Y(i)=c*U(i-1)+d*U(i-2)-a*Y(i-1)-b*Y(i-2);end; % creation de la sortie 

bruit_faible_amplitude= 0.1*randi([0 1],1,N); 

Y=Y+bruit_faible_amplitude;                             %Ajout de la bruit à la sortie
%-------------Question 3----------------------------------
% Identification par moindres carrés recursifs
%initialisation 
alpha=100;

%ou bien on peut initialiser teta par lea premier valeur 
A=[-Y(5) -Y(4) U(5) U(4);-Y(4) -Y(3) U(4) U(3);-Y(3) -Y(2) U(3) U(2);-Y(2) -Y(1) U(2) U(1)];
B=[Y(6);Y(5);Y(4);Y(3)];

teta=[10;17;20;-10]
%teta=pinv(A)*B;
P=alpha*eye(4);
TETA=teta;

%odre d'itiration = nombres des échantillons 
for k=7:N;   % il faut que K > na,nb+1
     h=[-Y(k-1) -Y(k-2) U(k-1) U(k-1)]';
     G    =         P*h*inv(1+h'*P*h);
     teta =         teta+G*(Y(k)-h'*teta);
     P    =         (eye(4)-G*h')*P;
     TETA =         [TETA,teta];
    
end;

display("\n les parametres estimés sont: \n")
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
    