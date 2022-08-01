"""
Objetivo: asignar los procedimientos
a la sala de tal forma que no hayan
cruces en los horarios seleccionados y se maximice
el tiempo que la sala se encuentre en
funcionamiento en el día.
"""

# Se abre y lee el archivo (entrada).
import numpy as np;
import math;
import time;

input = open("./sala_operaciones_entrada1.txt");
content = input.readlines();

"""
se asocian los indices de esta forma:
0 <- nombre del procedimiento.
1 <- hora de inicio del procedimiento.
2 <- hora de finalización del procedimiento.
content[2].split()[1]
"""

# funciones auxiliares
def merge(nam, tim, arr, l, m, r):
    n1 = m - l + 1
    n2 = r - m
 
    # create temp arrays, L1 -> names, L2 -> times (final or initial)
    L = [0] * (n1)
    R = [0] * (n2)
    L1 = [0] * (n1)
    R1 = [0] * (n2)
    L2 = [0] * (n1)
    R2 = [0] * (n2)
 
    # Copy data to temp arrays L[] and R[]
    for i in range(0, n1):
        L[i] = arr[l + i]
        L1[i] = nam[l + i]
        L2[i] = tim[l + i]
 
    for j in range(0, n2):
        R[j] = arr[m + 1 + j]
        R1[j] = nam[m + 1 + j]
        R2[j] = tim[m + 1 + j]
 
    # Merge the temp arrays back into arr[l..r]
    i = 0     # Initial index of first subarray
    j = 0     # Initial index of second subarray
    k = l     # Initial index of merged subarray

    while i < n1 and j < n2:
        if L[i] <= R[j]:
            nam[k] = L1[i]
            tim[k] = L2[i]
            arr[k] = L[i]
            i += 1
        else:
            nam[k] = R1[j]
            tim[k] = R2[j]
            arr[k] = R[j]
            j += 1
        k += 1
 
    # Copy the remaining elements of L[], if there
    # are any
    while i < n1:
        nam[k] = L1[i]
        tim[k] = L2[i]
        arr[k] = L[i]
        i += 1
        k += 1
 
    # Copy the remaining elements of R[], if there
    # are any
    while j < n2:
        nam[k] = R1[j]
        tim[k] = R2[j]
        arr[k] = R[j]
        j += 1
        k += 1
 
# l is for left index and r is right index of the
# sub-array of arr to be sorted
 
def mergeSort(nam, tim, arr, l, r):
    if l < r:
 
        # Same as (l+r)//2, but avoids overflow for
        # large l and h
        m = l+(r-l)//2
 
        # Sort first and second halves
        mergeSort(nam, tim, arr, l, m)
        mergeSort(nam, tim, arr, m+1, r)
        merge(nam, tim, arr, l, m, r)

def counting_sort(arr, max_value, get_index):
    counts = [0] * max_value

    # Counting - O(n)
    for a in arr:
        counts[get_index(a)] += 1;
  
    # Accumulating - O(k)
    for i, c in enumerate(counts):
        if i == 0:
            continue;
        else:
            counts[i] += counts[i-1]

    # Calculating start index - O(k)
    for i, c in enumerate(counts[:-1]):
        if i == 0:
            counts[i] = 0
        counts[i+1] = c

    ret = [None] * len(arr)
    # Sorting - O(n)
    for a in arr:
        index = counts[get_index(a)]
        ret[index] = a
        counts[get_index(a)] += 1
  
    return ret;

def paddingElementsLeft(e, n, val):
    """
    convierte al valor en uno de longitud n añadiendo elementos definidos:
    val <- valor (str/numeric)
    e <- elementos (str)
    n <- número de elementos (int); 
    """
    res = "";
    val = str(val);
    k = len(val);
    if k >= n:
        return val;
    for i in range(0,abs(k-n)):
        res += e;
    res += val;
    return res;

def minToHrs(m):
    if m > 1440:
        raise Exception("Error!, el resultado es mayor que 24hrs");
    hrs = math.trunc(m/60);
    min = (m - (hrs * 60));
    return "{hours}:{minutes}".format(hours = paddingElementsLeft("0",2,hrs), 
                                      minutes = paddingElementsLeft("0",2,min));

def hrsToMin(t):
    """
    formato de tiempo permitido: 00:00
    """
    h = int(t[0:2]);
    m = int(t[3:5]);
    res = h*60 + m;
    if res > 1440:
        raise Exception("Error!, el resultado es mayor que 24hrs");
    return res;

def timeSum(t1, t2):
    """
    formato de tiempo permitido: hh:mm
    """
    t1 = hrsToMin(t1);
    t2 = hrsToMin(t2);
    res = abs(t1 + t2);
    return minToHrs(res);

def timeDiff(t1, t2):
    """
    formato de tiempo permitido: hh:mm
    """
    t1 = hrsToMin(t1);
    t2 = hrsToMin(t2);
    res = abs(t1 - t2);
    return minToHrs(res);

# O(nlgn)
def sort_time(A, n):
    """
    ti -> arreglo con tiempos iniciales;
    tf -> arreglo con tiempos finales;
    n -> número de procedimientos;
    """
    # Hallamos el arreglo de los nombres de los procedimientos - O(n)
    names = [];
    for i in range(1, n+1):
        names.append(A[i].split()[0]);

    # Hallamos el arreglo de los tiempos iniciales (en minutos) - O(n)
    ti = [];
    for i in range(1, n+1):
        ti.append(hrsToMin(A[i].split()[1]));

    # Hallamos el arreglo de los tiempos finales (en minutos) - O(n)
    tf = [];
    for i in range(1, n+1):
        tf.append(hrsToMin(A[i].split()[2]));
    
    # Se ordenan los tiempos iniciales de menor a mayor - O(nlgn)
    # merge-sort tiene complejidad O(nlgn).
    mergeSort(names, tf, ti, 0, n-1);

    # Se convierten en horas de nuevo los valores de los procedimientos - O(n)
    for i in range(0, n):
        ti[i] = minToHrs(ti[i]);

    for i in range(0, n):
        tf[i] = minToHrs(tf[i]);

    # Se crea el arreglo final ordenado - O(n)
    sorted = [];
    for i in range(0, n):
        proc = names[i];
        t1 = ti[i];
        t2 = tf[i];
        sorted.append([proc,t1,t2]);

    return sorted;

# Se inicializan variables 
n = int(content[0].split()[0])

start = time.time();

# Cuerpo del algoritmo solución.
# Se ordenan los procedimientos siguiendo idea del radix-sort
sorted_content = sort_time(content, n);

# obtenemos la matriz del tiempo requerido para cada procedimiento (en minutos)
maxP = 0;
maxW = 0;
m = [];
for i in range(0,n):
    wp = hrsToMin( timeDiff(sorted_content[i][1], sorted_content[i][2]) );
    if maxW < wp:
        maxP = i;
        maxW = wp;
    m.append(wp);
# print("m:", m,"maxP",maxP,"maxW",maxW);

# Se comparan los tiempos finales del procedimiento actual con el tiempo inicial 
# del siguiente para determinar la solución voraz.

# obtengo una solución partiendo desde donde esta el inicio, con su valor.
sol1 = [0];
w1 = m[0];
k1 = 0;
for i in range(1,n):
    cw = m[i];
    pw = m[sol1[k1]];
    p2w = m[sol1[k1-1]];
    cit = hrsToMin(sorted_content[i][1]);
    pft = hrsToMin(sorted_content[sol1[k1]][2]);
    p2ft = hrsToMin(sorted_content[sol1[k1-1]][2]);
    if pft <= cit:
        k1 += 1;
        w1 += cw;
        sol1.append(i);
    else:
        if p2ft <= cit and cw > pw:
            w1 += - pw + cw;
            sol1[k1] = i;

# obtengo una solución desde donde esté el procedimiento con peso máximo, con su valor.
sol2 = [maxP]
w2 = maxW;
k2 = 0;
for i in range(maxP+1,n):
    cw = m[i];
    pw = m[sol2[k2]];
    p2w = m[sol2[k2-1]];
    cit = hrsToMin(sorted_content[i][1]);
    pft = hrsToMin(sorted_content[sol2[k2]][2]);
    p2ft = hrsToMin(sorted_content[sol2[k2-1]][2]);
    if pft <= cit:
        k2 += 1;
        w2 += cw;
        sol2.append(i);
    else:
        if p2ft <= cit and cw > pw:
            w2 += - pw + cw;
            sol2[k2] = i;
# print("sol1",sol1,"w1",w1,"sol2",sol2,"w2",w2)

sol = [];
timeMax = 0;
c = 0;
if w2 > w1:
    sol = sol2;
    timeMax = w2
    c = k2+1;
else:
    sol = sol1;
    timeMax = w1;
    c = k1+1;
# print("sol",sol,"c",c,"timeMax",timeMax)

end = time.time();

print("--- execution time: {t} ---".format(t=end-start));

input.close();

# Se genera la salida en formato .txt
output = open("./sala_operaciones_salida.txt", "w");
output.write(f"{c} -- procedimientos\n");
output.write("{time} -- tiempo de uso\n".format(
    time=timeMax));
for i in range(0,c):
    output.write(f"{sorted_content[sol[i]][0]}\n");
output.close();