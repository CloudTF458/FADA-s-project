"""
Objetivo: Su trabajo consiste en determinar como repartir
los libros entre los escritores de forma que
el tiempo que se demoren realizando las copias
sea el mínimo teniendo en cuenta que el tiempo
que se demoran en realizar las copias depende
del escritor que tenga mas paginas asignadas.
Hacer la copia de una pagina tiene una duracion de un dia.

NOTA: un escritor le puede
corresponder los libros que estan entre l0 y l3 a
otro escritor desde l4 y l6 y así sucesivamente.
"""

# se importa libreria math
import numpy as np;
import time;
import copy;
import math;

# Se abre y lee el archivo (entrada).
input = open("./libros_entrada2.txt");
content = input.readlines();

"""
se asocian los indices de esta forma:
0 <- nombre del libro.
1 <- numero de páginas del libro.
content[2].split()[1]
"""

# funciones auxiliares
def translate(W, sol):
    """
    O(n)
    devuelve el arreglo de pesos asociados a los indices del arreglo solución
    ejemplo:
    entrada: [1, 4, 5]
    salida: [275.0, 335.0, 350.0]
    Esto significa que para un total de 6 libros (considerando que el máximo indice es 5):
    el autor 1 toma 2 libros que tienen 275 páginas,
    el autor 2 toma 3 libros que tienen 335 páginas y 
    el autor 3 toma 1 libros que tiene 350 páginas.
    """
    sols = [];
    sols.append(W[0][sol[0]]);
    for i in range(1, n):
        c = abs(sol[i-1] - sol[i]) - 1;
        sols.append(W[sol[i]-c][sol[i]]);
    return sols;

def voraz_sol(W, n, m):
    """
    entradas:
    W -> matriz con número de paginas de cada libro.
    n -> número de escritores.
    m -> número de libros.
    
    devuelve el rango de libros que se asignará a cada escritor en la 
    solución optima.
    ej: [3, 4, 6]
    significa:
    autor 1 -> lib 1 - lib 3
    autor 2 -> lib 4 - lib 5
    autor 3 -> lib 6
    """
    
    prom = 0;
    for i in range(0, m):
        prom += W[i];
    prom = prom/m

    sd = 0;
    for i in range(0, m):
        sd += ((W[i]-prom)**2);
    sd = (sd/m)**(1/2);

    print("prom",prom,"sd",sd)

    sol = [];
    if sd == 0:
        c = m-1;
        print("c",c)
        for i in range(1, n+1):
            maxI = math.ceil(c/n);
            sol.append(c);
            c -= maxI;
        fSol = sol;
        sol = [];
        for i in range(0, n):
            sol.append(fSol[n-i-1]);
    else:
        acum = 0;
        for i in range(0, m):
            print("acum",acum,"i",i,"sol",sol)
            acum += W[i];
            if acum > maxP:
                acum = 0;
                sol.append(i);
            if len(sol) == n-1:
                break;
        sol.append(m-1);

    return sol
            

# Se inicializan variables 
# Número de escritores.
n = int(content[0].split()[0])
# Número de libros.
m = int(content[0].split()[1]);

# Número de libros real
nBooks = len(content) - 1;

start = time.time();

# Cuerpo del algoritmo solución.
if not n > m and not n < 1 and m == nBooks:
    # Se crea vector de paginas para usarlo en la creación de la matriz de sumas
    pages = [];
    maxL = 0;
    maxP = int(content[1].split()[1]);
    pages.append(maxP);
    for i in range(2,m+1):
        pag = int(content[i].split()[1])
        if maxP < pag:
            maxP = pag;
            maxL = i-1;
        pages.append(pag);
    print("pages",pages);

    # Se crea matriz de sumas
    sums = np.zeros((m,m));
    c = 0;
    for i in range(0,m):
        c = 0;
        for j in range(i,m):
            c += pages[j];
            sums[i][j] = c;

    # print("pages:",pages)
    print("sums:")
    print(sums);

    # Se halla la solución optima para los lectores.
    sol = voraz_sol(pages, n, m);
    maxTime = max(translate(sums, sol));

    print("solFinal:",sol,"t",translate(sums,sol))
    print("pags")
    print(maxTime)
else:
    raise Exception("Error! invalid input values:",(n,m));

end = time.time();

print(" --- execution time: {t} --- ".format(t = end-start));

input.close();

# Se genera la salida en formato .txt
# output = open("./libros_salida.txt", "w");
# output.write(f"{maxTime} -- dias requeridos\n");
# if sol[0] > 0:
#     output.write(f"lib{1} - lib{sol[0]+1}\n");
# else:
#     output.write(f"lib{sol[0]+1}\n");
# if n > 1:
#     for i in range(1,n-1):
#         # print("i",i)
#         if (sol[i] - sol[i-1]) != 1:
#             output.write(f"lib{sol[i-1]+2} - lib{sol[i]+1}\n");
#         else:
#             output.write(f"lib{sol[i]+1}\n");
#     if (sol[n-1] - sol[n-2]) != 1:
#         output.write(f"lib{sol[n-2]+2} - lib{m}\n");
#     else:
#         output.write(f"lib{sol[n-1]+1}");
# output.close();