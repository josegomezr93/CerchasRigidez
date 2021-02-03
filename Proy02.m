%% Datos de entradas
%Script que sirve para el calculo de cerchas ingresando las coordenadas de
%los nodos, propiedades de los materiales, restricciones en el id del nodo
%correspondiente y cargas externas

clc, clear all, close all
nNodos = input('Introducir el numero de nodos que conforman la cercha: '); %Numero de nodos que conforman la cercha
nElementos = input('Numero de elementos: '); %Numero de elementos

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
Es = zeros(nElementos,2);
As = zeros(nElementos,2);
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
    Es(i,1) = i;
    As(i,1) = i;
    Es(i,2) = input('Introducir Modulo de elasticidad del elemento: ');
    As(i,2) = input('Introducir Area de la seccion: ');
    longEle(i,2) = sqrt(power(dX,2)+power(dY,2));
    cosEle(i) = dX./longEle(i,2);
    senEle(i) = dY./longEle(i,2);
end

%% Calculo de Grados de Libertad
idR = input('Introducir el Id de las coordenadas globales que se restriguen: ');
Rest = length(idR); %Numero de restricciones en nodos
GDL = (2*nNodos) - Rest;
display(GDL);

%% Matrices de rigidez de cada elemento
k = zeros(2,2, nElementos);
eleContador = 1;
while eleContador <= nElementos;
    i = 1;
    while i < 2
        j = i;
        k(i,j,eleContador) = Es(eleContador,2)*As(eleContador,2)/longEle(eleContador,2);
        k(i, j+1,eleContador) = -Es(eleContador,2)*As(eleContador,2)/longEle(eleContador,2);
        k(i+1,j+1,eleContador) = k(i,j,eleContador);
        k(i+1,j,eleContador) = k(i, j+1,eleContador);
        i = i+1;
    end
    eleContador = eleContador +1;
end
%%
%Matriz de transformacion T
a = 4; b = 2; %Variables que definen las dimensiones del vector Transformacion
T = zeros(a,b,nElementos); %Matriz de transformacion para guardarla en un solo espacio de memoria
for i = 1:nElementos
    T(1,1,i) = cosEle(i);
    T(2,1,i) = senEle(i);
    T(3,2,i) = T(1,1,i);
    T(4,2,i) = T(2,1,i);
end

%Matriz de Rigidez Transformada a coordenadas Globales
aux = 4;
kG = zeros(aux,aux,nElementos);
for i = 1:nElementos
    kG(:,:,i) = T(:,:,i)*(k(:,:,i)*T(:,:,i)');
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

%%
aux = 1;
gdlNodo = 4; %GDL por barra || Se consideran dos coordenadas por nodo
dimXY = gdlNodo + 1; %Se adiciona una fila y columna que contendra el id de la coordenada global
GlobalK = zeros(dimXY,dimXY,nElementos);
coordElemIJ = zeros(1,dimXY);
uEle = zeros(2,2,nElementos); %Instancio la variable donde se va guardaran los resultados
pEle = uEle;
iMatriz = 1;
while aux <= nElementos
    auxI = [2*idEleNodos(aux,2)-1 2*idEleNodos(aux,2)];
    auxJ = [2*idEleNodos(aux,3)-1 2*idEleNodos(aux,3)];
    GlobalK(1,:,aux) = [0 auxI auxJ];
    GlobalK(:,1,aux) = GlobalK(1,:,aux)';
    GlobalK(2:end,2:end,aux) = kG(:,:,aux);
    auxNodo = [idEleNodos(aux,2) idEleNodos(aux, 3)];
    uEle(:,1,aux) = auxNodo';
    pEle(:,1,aux) = auxNodo';
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

Ue = zeros(gdlNodo,2,nElementos); %Desplazamientos nodales Coordenadas globales
peGlobal = zeros(gdlNodo, 2, nElementos); %Esfuerzos en las barras en Coordenadas globales
i = 1; %Variable auxiliar que sirve de contador
iMatriz = 1; %Variable auxiliar que sirve para armar los vectores de esfuerzos en las barras

while i <= nElementos
        auxI = [2*idEleNodos(i,2)-1 2*idEleNodos(i,2)];
        auxJ = [2*idEleNodos(i,3)-1 2*idEleNodos(i,3)];
		idUe = [auxI auxJ];
        peGlobal(:,1,i) = idUe';
        Ue(:,1,i) = idUe';
        Ue(:,2,i) = D([Ue(:,1,i)],2);
        peGlobal(:,2,i) = kG(:,:,i)*Ue(:,2,i);
        uEle(:,2, i) = T(:,:,i)'*Ue(:,2,i); %Deformaciones en las barras Coordenadas Locales
        pEle(:,2,i) = k(:,:,i)*uEle(:,2,i); %Calculo de esfuerzo en las barras Coordenadas Locales
        i = i + 1;
end

%%
%Reacciones 
pR = zeros(2*nNodos, 2); %Vector donde se guardaran las reacciones en los nodos correspondientes
fR = 0; %Variable auxiliar para guardar en cada paso de tiempo el valor de las fuerzas internas dentro de cada elemento
pR(:,1) = 1:2*nNodos;
j = 1;
while j <= nElementos
    for r = 1:length(peGlobal(:,1,1))
            auxR = peGlobal(r,1,j); %Ubicacion de la posicion de la hipermatriz con respecto a las fuerzas internas de cada elemento
            fR = peGlobal(r,2,j); %Actualizo vector de fuerza interna para utilizar y calcular las reacciones
            pR(auxR,2) = pR(auxR,2) + fR;
    end
        j = j + 1;
end
pR = pR(idR',:); %Actualizo vector de reacciones de acuerdo a las coordenadas restringidas

%% Presentacion de los resultados
%Tabla resumen resultados esfuerzos internos y deformaciones en cada barra
resDef = zeros(nElementos,3); %Matriz resumen
resEsf = zeros(nElementos, 5);
resDef(:,1) = 1:nElementos;
resEsf(:,1) = resDef(:,1);
for i = 1:nElementos
    resEsf(i,2) = pEle(1,2,i);
    resEsf(i,3) = pEle(2,2,i);
    resDef(i,2) = Ue(1,2,i);
    resDef(i,3) = Ue(2,2,i);
    resDef(i,4) = Ue(3,2,i);
    resDef(i,5) = Ue(4,2,i);
end
 
tablaEsf = table([resEsf(:,1), resEsf(:,2), resEsf(:,3)],'variableNames', {['idEle ', 'N_i [Kgf]', ' N_j [Kgf]']});
tablaDef = table([resDef(:,1), resDef(:,2), resDef(:,3), resDef(:,4), resDef(:,5)],'variableNames', {['idEle ', 'U_ix [m]', ' U_iy [m] ', 'U_jx [m]', ' U_jy [m]']});
tablaReacciones = table([pR(:,1), pR(:,2)], 'variableNames', {['GDL ', ' Reac [Kgf]']});
display(tablaEsf); display(tablaDef); display(tablaReacciones);
