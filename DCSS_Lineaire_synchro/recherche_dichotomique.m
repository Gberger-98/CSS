function [pos] = recherche_dichotomique(a, b, epsilon, x)

N = numel(x);

H = @(w)(abs( sum(x.*exp(-1i.*(0:N-1)*w)) ));

p = log2((b-a)/epsilon) + 1;
M = 2^p;

x_gauche = a;
x_droite = b;
pas = (b-a)/M;
y_gauche = H(a);
y_droite = H(b);

for i = 1:p
    y_g = H((x_droite+x_gauche)/2);
    y_d = H((x_droite+x_gauche)/2 + pas);
    A = [y_gauche, y_g, y_d, y_droite];
    pos = find(A==max(A));
    if((pos(1)==1) || (pos(1)==2))
        y_droite = y_g;
        x_droite = (x_droite+x_gauche)/2;
    else
        y_gauche = y_d;
        x_gauche = (x_droite+x_gauche)/2+pas;
    end
end

A = [H(x_gauche),H(x_droite)];
temp = find(A==max(A));
if (temp==1)
    pos = x_gauche;
else
    pos = x_droite;
end


