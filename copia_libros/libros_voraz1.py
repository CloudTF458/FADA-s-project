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
import copy
import time;

# Se abre y lee el archivo (entrada).
input = open("./libros_entrada8.txt");
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

def optimal_sol(W, sol, c, maximum, best, prevsol):
    """
    devuelve el rango de libros que se asignará a cada escritor en la 
    solución optima.
    ej: [3, 4, 6]
    significa:
    autor 1 -> lib 1 - lib 3
    autor 2 -> lib 4 - lib 5
    autor 3 -> lib 6
    """
    sols = translate(W, sol);
    best = max(sols);
    # print("c",c,"n",n,"best",best,"max",maximum,"sol",sol,"sols",sols,"prevsol",prevsol)
    if c==n:
        return sol;
    if best < maximum:
        maximum = best;
        prevsol = sol.copy();        
        c = 0;
    for i in range(0,n):
        if sols[i] == best:
            # print("sol[i]",sol[i],"sol[i-1]",sol[i-1],"w:",sols[i]);
            c += 1;
            if i == 0 and sol[i] == 0:
                return sol;
            if i == n-1:
                return prevsol;
            if (sol[i] - sol[i-1]) == 1:
                return sol;
            sol[i] = sol[i] - 1;
            # print("best",best,"max",maximum)
            return optimal_sol(W, sol, c, maximum, best, prevsol);

def voraz_sol(W, p, n, m):
    """
    entradas:
    W -> matriz de sumas.
    p -> matriz con número de paginas de cada libro.
    n -> número de escritores.
    m -> número de libros.

    O(n*m)
    devuelve el rango de libros que se asignará a cada escritor en la 
    solución optima.
    ej: [3, 4, 6]
    significa:
    autor 1 -> lib 1 - lib 3
    autor 2 -> lib 4 - lib 5
    autor 3 -> lib 6
    """

    iSol = [];
    wsol = 0;

    # Se llena una solución, cuyo ultimo indice es m-1 - O(n)
    
    if True:
        for i in range(m-n, m):
            iSol.append(i);
    else:
        for i in range(0, n-1):
            iSol.append(i);
        iSol.append(m-1);

    best = max(translate(W,iSol))
    bSol = copy.deepcopy(iSol);
    sol = copy.deepcopy(iSol);

    print("bsol",bSol,"w",best);

    # Se busca la solución voraz recorriendo la matriz de sumas
    
    for i in range(0, n-1):
            if i == 0:
                for j in range(i, bSol[i]):
                    sol[i] = j;
                    wsol = max(translate(W,sol))
                    print("i",i,"j",j,"sol:",sol,"w",translate(sums,sol))
                    if wsol < best:
                        best = wsol;
                        bSol = sol;
            else:
                for j in range(bSol[i-1], bSol[i+1]):
                    sol[i] = j;
                    wsol = max(translate(W,sol))
                    print("i",i,"j",j,"sol:",sol,"w",translate(sums,sol))
                    if wsol < best:
                        best = wsol;
                        bSol = sol;
    return bSol

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
    for i in range(1,m+1):
        pages.append(int(content[i].split()[1]));
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

    # Se halla la mejor solución (optima) para los lectores.
    first_sol = [];
    for i in range(0,n):
        first_sol.append(m-n+i);

    first_sol = np.array(first_sol);
    sol = optimal_sol(sums, first_sol, 0, max(translate(sums, first_sol)), 0, first_sol.copy());
    maxTime = max(translate(sums, sol));

    print("sol:")
    print(sol);
    print("pags")
    print(translate(sums, sol))
else:
    raise Exception("Error! invalid input values:",(n,m));

end = time.time();

print(" --- execution time: {t} --- ".format(t = end-start));

input.close();

# Se genera la salida en formato .txt
output = open("./libros_salida.txt", "w");
output.write(f"{maxTime} -- dias requeridos\n");
if sol[0] > 0:
    output.write(f"lib{1} - lib{sol[0]+1}\n");
else:
    output.write(f"lib{sol[0]+1}\n");
if n > 1:
    for i in range(1,n-1):
        # print("i",i)
        if (sol[i] - sol[i-1]) != 1:
            output.write(f"lib{sol[i-1]+2} - lib{sol[i]+1}\n");
        else:
            output.write(f"lib{sol[i]+1}\n");
    if (sol[n-1] - sol[n-2]) != 1:
        output.write(f"lib{sol[n-2]+2} - lib{m}\n");
    else:
        output.write(f"lib{sol[n-1]+1}");
output.close();