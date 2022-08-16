%Partie_2
%Initialisation 

display ("\nles parametres sont :\n")
a=1
b=0.7
c=0.8
d=-1
e=1.5
f=1



N=100 ;

U=randi([0 1],1,N); 

% Y=zeros(1,N);
Y(1)=U(1);
Y(2)=U(2);
Y(3)=U(3);
%cas où les paramètres a, b, c et d ne sont pas constants 
% mais convergents vers des valeurs constantes.
for i=4:N;
    a_k= a*(1-0.5*exp(-i));
    b_k= b*(1-0.5*exp(-i));
    c_k= c*(1-0.5*exp(-i));
    
    d_k= d*(1-0.5*exp(-i));
    e_k= e*(1-0.5*exp(-i));
    f_k= f*(1-0.5*exp(-i));



    Y(i)=d_k*U(i-1)+e_k*U(i-2)+f_k*U(i-3)-a_k*Y(i-1)-b_k*Y(i-2)-c*U(i-3);
end; % creation de la sortie 

bruit_faible_amplitude= 0.1*randi([0 1],1,N); 



display ("   1ère cas : avec bruit \n")

%----------------------------------------------------------


% Identification par moindres carrés recursifs
%initialisation 
alpha=1000;

%ou bien on peut initialiser teta par lea premier valeur 
A=[-Y(8) -Y(7) -Y(6) U(8) U(7) U(6)];

for i=1:5;
l=[-Y(8-i) -Y(7-i) -Y(6-i) U(8-i) U(7-i) U(6-i)];
A=[A;l];
end;

B=[Y(8);Y(7); Y(6); Y(5); Y(4); Y(3)];

teta=pinv(A)*B;
P=alpha*eye(6);
TETA1=teta;
% pour le cas ou a,b,c et d sont variables on utilise l'algorithme avec
% facteur d'oublie
%odre d'itiration = nombres des échantillons 


lambda = .95 ; % facteur d'oublie 

for k=10:N;  
     h=[-Y(k-1); -Y(k-2) ;-Y(k-3); U(k-1); U(k-2) ;U(k-3)];
     G    =         P*h*inv(lambda+h'*P*h);
     teta =         teta+G*(Y(k)-h'*teta);
     P    =         (1/lambda)*(eye(6)-G*h')*P;
     TETA1 =         [TETA1,teta];
end;
a_nb=teta(1)
b_nb=teta(2)
c_nb=teta(3)

d_nb=teta(4)
e_nb=teta(5)
f_nb=teta(6)

%--------------------------------------------------
display ("\n    2ème cas : avec bruit\n") 
Y=Y+bruit_faible_amplitude;


% Identification par moindres carrés recursifs
%initialisation 
alpha=1000;

%ou bien on peut initialiser teta par lea premier valeur 
A=[-Y(8) -Y(7) -Y(6) U(8) U(7) U(6)];

for i=1:5;
l=[-Y(8-i) -Y(7-i) -Y(6-i) U(8-i) U(7-i) U(6-i)];
A=[A;l];

end;

B=[Y(8);Y(7); Y(6); Y(5); Y(4); Y(3)];

teta=pinv(A)*B;
P=alpha*eye(6);
TETA=teta;
% pour le cas ou a,b,c et d sont variables on utilise l'algorithme avec
% facteur d'oublie
%odre d'itiration = nombres des échantillons 


lambda = .95;  % facteur d'oublie 

for k=10:N;  
     h=[-Y(k-1); -Y(k-2) ;-Y(k-3); U(k-1); U(k-2) ;U(k-3)];
     G    =         P*h*inv(lambda+h'*P*h);
     teta =         teta+G*(Y(k)-h'*teta);
     P    =         (1/lambda)*(eye(6)-G*h')*P;
     TETA =         [TETA,teta];
end;


a_b=teta(1)
b_b=teta(2)
c_b=teta(3)
d_b=teta(4)
e_b=teta(5)
f_b=teta(6)
% calcule de la marge d'erreur pour les deux cas 

Erreur_san_bruit=[a;b;c;d;e;f]-[a_nb;b_nb;c_nb;d_nb;e_nb;f_nb]

Erreur_avec_bruit=[a;b;c;d;e;f]-teta
