%% Datos de entradas
clc, clear all, close all
nNodos = 4; %Numero de nodos que conforman la cercha
% nElementosPrincipales = 31;
% nElementosSecundarios = 18;
% nElementosTotales = nElementosPrincipales + nElementosSecundarios; %Numero de elementos totales que conforman la cercha
nElementos = 6; %Numero de elementos
% ACordonesBorde = 7610; %Area para los cordones de borde [cm2]
% AZonaApoyo = 9480; %Area para las barra que conforman el triangulo en la zona de apoyo [cm2]
% ASecundarias = 900; %Area para las barras del tejido [cm2]
As = 100; %area seccion transversal de los elementos [cm2]

Es = 2500; %Modulo de elasticidad del acero [tonf/cm2]
%Em = 1.5e5; %Modulo de elasticidad de la madera (tejido)

%Coordenadas de los nodos
coordNodos = zeros(nNodos,3);
for i = 1:nNodos
    fprintf('Coordenadas del nodo %i \n',i);
    coordNodos(i,1) = i; %Identificador del nodo
    coordNodos(i,2) = input('Introducir dX del nodo: ');
    coordNodos(i,3) = input('Introducir dY del nodo: ');
end

idEleNodos = zeros(nElementos,3);
longEle = zeros(nElementos,2);
cosEle = zeros(1, nElementos);
senEle = zeros(1, nElementos);
for i = 1:nElementos
    fprintf('Elemento %i \n', i);
    nodoI = input('Introducir el id del nodo inicial del elemento: ');
    nodoJ = input('Introducir el id del nodo final del elemento: ');
    idEleNodos(i,1) = i; %id del elemento
    idEleNodos(i,2) = nodoI; %id del nodo I para el elemento
    idEleNodos(i,3) = nodoJ; %id del nodo J para el elemento
    dX = coordNodos(nodoJ,2) - coordNodos(nodoI,2);
    dY = coordNodos(nodoJ,3) - coordNodos(nodoI,3);
    longEle(i,1) = i; 
    longEle(i,2) = sqrt(power(dX,2)+power(dY,2));
    cosEle(i) = dX./longEle(i,2);
    senEle(i) = dY./longEle(i,2);
end

display(longEle);  display(cosEle); display(senEle);

%% Calculo de Grados de Hipergeometria
idR = input('Introducir el Id de las coordenadas globales que se restriguen: ');
Rest = length(idR); %Numero de restricciones en nodos
GDL = (2*nNodos) - Rest;
display(GDL);

%% Matrices de rigidez de cada elemento
k = zeros(2*nElementos);
ele = 1
for i = 1:2:(2*nElementos)
    j = i;
    k(i,j) = Es*As/longEle(ele,2);
    k(i, j+1) = -Es*As/longEle(ele,2);
    k(i+1,j+1) = k(i,j);
    k(i+1,j) = k(i, j+1);
    ele = ele +1;
end

%Matriz de transformacion T
a = nElementos*4; b = nElementos*2;
T = zeros(a,b); %Matriz de transformacion para guardarla en un solo espacio de memoria
ele = 1;
j = 1
for i = 1:4:a
    T(i,j) = cosEle(ele);
    T(i+1,j) = senEle(ele);
    T(i+2,j+1) = T(i,j);
    T(i+3,j+1) = T(i+1,j);
    j = j + 2;
    ele = ele + 1;
end

%Matriz de Rigidez Transformada a coordenadas Globales
aux = 4*nElementos;
kG = zeros(aux);
aux1 = 1;
aux2 = 1;
auy1 = 1;
for i = 1:4:aux
    kG(i:3+i,i:3+i) = T(aux1:aux1+3,auy1:auy1+1)*(k(aux2:aux2+1,aux2:aux2+1)*T(aux1:aux1+3,auy1:auy1+1)')
    aux1 = aux1 + 4;
    auy1 = auy1 + 2;
    aux2 = aux2 + 2;
end

%Vector de restricciones R
R = zeros(2*nNodos,2);
F = zeros(2*nNodos,2);
for i = 1:2*nNodos
    R(i,1) = i; %Identificador de la coordenada que se esta evaluando
    F(i,1) = i;
    if find(i == idR)
        R(i,2) = 1; %Si se restringue esa coordenada global
    else
        R(i,2) = 0; %No se restringue esa coordenada global
    end
end

aux = 1;
gdlNodo = 4; %GDL por barra || Se consideran dos coordenadas por nodo
dimXY = gdlNodo + 1; %Se adiciona una fila y columna que contendra el id de la coordenada global
GlobalK = zeros(dimXY,dimXY,nElementos);
coordElemIJ = zeros(1,dimXY);
iMatriz = 1;
while aux <= nElementos
    auxI = [2*idEleNodos(aux,2)-1 2*idEleNodos(aux,2)];
    auxJ = [2*idEleNodos(aux,3)-1 2*idEleNodos(aux,3)];
    GlobalK(1,:,aux) = [0 auxI auxJ];
    GlobalK(:,1,aux) = GlobalK(1,:,aux)';
    GlobalK(2:end,2:end,aux) = kG(iMatriz:iMatriz+3,iMatriz:iMatriz+3);
    iMatriz = iMatriz + 4;
    aux = aux +1;
end

%Ensamblaje de la matriz de rigidez K
K = zeros(GDL);
R([find(R(:,2)==1)],:) = [];
R(:,2) = [];
for i = 1:GDL
    aux1 = R(i); %Indice de fila
    for j = 1:GDL
        aux2 = R(j); %Indice de columna
        m = 1;
        while m <= nElementos
            if find(aux1 == GlobalK(2:end,1,m) & aux2 == GlobalK(1,2:end,m))
                [indeX indeY] = find(aux1 == GlobalK(2:end,1,m) & aux2 == GlobalK(1,2:end,m));
                K(i,j) = K(i,j) + GlobalK(indeX+1,indeY+1,m);
            end
            m = m+1;
        end
    end
end
            
%% Vector Cargas
for i = 1:2*nNodos
    fprintf('Coordenada Global %i\n', i);
    F(i,2) = input('Introducir valor de carga externa correspondiente a la coordenada Global: ');
end
for i = 1:length(idR)
    aux1 = find(F(:,1) == idR(i));
    F(aux1,:)=[];
end
d = K\F(:,2); %Vector Desplazamientos Globales para los GDL
dGlobal = zeros(length(d),2);
dGlobal(:,1) = R;
dGlobal(:,2) = d;
%%
D = zeros(2*nNodos, 2);
for i = 1:2*nNodos
    D(i,1) = i;
    if find(D(i,1) == dGlobal(:,1))
        id = find(D(i,1) == dGlobal(:,1));
        D(i,2) = dGlobal(id,2);
    end
end

Ue = zeros(gdlNodo,2,nElementos);
peGlobal = zeros(gdlNodo, 2, nElementos);
i = 1;
iMatriz = 1;
while i <= nElementos
        auxI = [2*idEleNodos(i,2)-1 2*idEleNodos(i,2)];
        auxJ = [2*idEleNodos(i,3)-1 2*idEleNodos(i,3)];
		idUe = [auxI auxJ];
        peGlobal(:,1,i) = idUe';
        Ue(:,1,i) = idUe';
        Ue(:,2,i) = D([Ue(:,1,i)],2);
        peGlobal(:,2,i) = kG(iMatriz:iMatriz+3,iMatriz:iMatriz+3)*Ue(:,2,i);
        iMatriz = iMatriz + 4;
        i = i + 1;
end
%%
%Reacciones 
pR = zeros(length(idR), 2);
fR = zeros(2*nNodos,2);
fR(:,1) = 1:2*nNodos;
pR(:,1) = idR';
for i = 1:length(fR)
    auxR = idR(i);
    if find(fR(i,1) == auxR)
        aux = find(fR(i,1) == auxR);
        fR(i,2) = F(auxR, 2);
    end
end

for i = 1:length(idR)
    auxR = idR(i);
    j = 1;
    while j <= nElementos
        if find(peGlobal(:,1, j) == auxR)
            aux = find(peGlobal(:,1, j) == auxR);
            pR(i,2) = peGlobal(aux,2,j) + pR(i,2); %Calculo el vector de reacciones
        end
        j = j + 1;
    end
    %fextR = find(fR(:,1) == auxR);
    %pR(i,2) = pR(i,2) + fextR;
end

%% Esfuerzos en los elementos


