import pydicom as dc
import numpy as np

#Extraer los datos
data = dc.dcmread("DET2EXTR001_DS.dcm", force=True)
data = data.pixel_array
#data = np.array(((0,0,0,0,0,0,0,0,0,0,0,0,0,0),(0,0,0,0,0,0,0,0,0,0,0,0,0,0),(0,0,1,1,1,1,1,9,8,1,1,1,1,1),(0,0,2,2,2,2,2,2,4,4,4,4,4,4),(0,0,3,3,3,3,3,3,4,4,4,4,3,3), (0,0,2,2,2,2,2,2,2,2,2,2,2,2), (0,0,2,2,2,2,2,0,5,5,5,5,5,5), (0,0,2,2,4,4,4,4,4,4,2,2,2,2), (0,0,1,1,1,1,1,9,8,1,1,1,1,1),(0,0,2,2,2,2,2,2,6,0,6,6,6,6), (0,0,1,1,1,1,1,3,3,1,1,1,1,1),(0,0,2,2,2,2,2,2,2,2,2,2,2,2),(0,0,0,0,0,0,0,0,0,0,0,0,0,0),(0,0,0,0,0,0,0,0,0,0,0,0,0,0),(0,0,0,0,0,0,0,0,0,0,0,0,0,0)))
print("Raw data")
print(data)
print("")

# restar filas de 0s
def borrar0(dicom):

    l = list(dicom.shape)
    filas = l[0]
    columnas = l[1]
    borrarf = list()
    borrarc = list()

    for x in range(0, filas):
        lista = dicom[x,:]
        if sum(lista) < 1:
            borrarf.append(x)

    for y in range(0, columnas):
        lista = dicom[:,y]
        if sum(lista) < 1:
            borrarc.append(y)
    
    d_sin0 = np.delete(dicom, borrarf, 0)
    d_sin0 = np.delete(d_sin0, borrarc, 1)

    return d_sin0
d1 = borrar0(data)
print("restar filas de 0s")
print(d1)
print("")

# funciones de suavizado y calculo de la uniformidad
def suavizado(d, f, c):
 
 # Restar bordes 
                              
 rf = range(0, f)
 lf = list(rf)
 rc = range(0, c)
 lc = list(rc)
 rf2 = range(filas-f, filas)
 lf2 = list(rf2)
 rc2 = range(columnas-c, columnas)
 lc2 = list(rc2)
                    
 m = np.delete(d, lf2, 0)
 m = np.delete(m, lc2, 1)
 m = np.delete(m, lf, 0)
 m = np.delete(m, lc, 1) 
 dfilt = np.copy(m)

 # Filtro

 lm = list(m.shape)
 filasm = lm[0]-1
 columnasm = lm[1]-1
 

 for x in range(0,filasm+1):
  for y in range(0, columnasm+1):
     if x == 0  or y == 0:
         a1 = 0
         b1 = 0
     else:
         a1 = m[x-1][y-1]
         b1 = 1
     if x == 0:
         a2 = 0
         b2 = 0
     else:
         a2 = 2*m[x-1][y]
         a2 =2
     if x == 0 or y == columnasm:
         a3 = 0
         b3 = 0
     else:   
         a3 = m[x-1][y+1]
         b3 = 1
     if y == 0:
         a4 = 0
         b4 = 0
     else:
         a4 = 2*m[x][y-1]
         b4 = 2

     a5 = 4*m[x][y]
     b5 = 4

     if y == columnasm:
         a6 = 0
         b6 = 0
     else:
         a6 = 2*m[x][y+1]
         b6 = 2

     if x == filasm or y == 0:
         a7 = 0   
         b7 = 0 
     else:
         a7 = m[x+1][y-1]
         a7 = 1

     if x == filasm:
         a8 = 0
         b8 = 0
     else:
         a8 = 2*m[x+1][y]
         b8 = 2
     if x == filasm or y == columnasm:
         a9 = 0
         b9 = 0
     else:
         a9 = (m[x+1][y+1])
         b9 = 1

     
     
     dfilt[x,y] = float((a1 + a2 + a3 + a4 + a5 + a6 + a7 + a8 + a9)/(b1 + b2 + b3 + b4 + b5 + b6 + b7 + b8 + b9))
     
 return dfilt
def uniformidad(matriz):
 
 # UFOV integral
 
 datosUFOV = matriz.flatten()
 listUFOV = list()
 for w in range (0, len(datosUFOV)):
    if datosUFOV[w] != 0:
       listUFOV.append(float(datosUFOV[w]))
 

 maximo = max(listUFOV)
 print(f"Cuentas maximas {maximo}")
 minimo = min(listUFOV)
 print(f"Cuentas minimas {minimo}")
 UFOVi = 100*(maximo - minimo)/(maximo + minimo)
 print(f"\n UFOV integral {UFOVi}")

 # CFOV integral

 l = list(matriz.shape)
 filas = l[0]
 columnas = l[1]
 restof=round(0.25*filas)
 restoc=round(0.25*columnas)
 
 rf = range(0, restof//2)
 lf = list(rf)
 rc = range(0, restoc//2)
 lc = list(rc)
 rf2 = range(filas-restof//2, filas)
 lf2 = list(rf2)
 rc2 = range(columnas-restoc//2, columnas)
 lc2 = list(rc2)
 
 mCFOV = np.delete(matriz, lf2, 0)
 mCFOV = np.delete(mCFOV, lc2, 1)
 mCFOV = np.delete(mCFOV, lf, 0)
 mCFOV = np.delete(mCFOV, lc, 1)
 
 datosCFOV = mCFOV.flatten()
 listCFOV = list()
 for w in range (0, len(datosCFOV)):
    if datosCFOV[w] != 0:
       listCFOV.append(float(datosCFOV[w]))

 maximo = max(listCFOV)
 minimo = min(listCFOV)
 CFOVi = 100*(maximo - minimo)/(maximo + minimo)
 print(f"\n CFOV integral {CFOVi}")

 # UFOV diferencial
   
 dif = 0
 for y in range(0,len(listUFOV),5):
      
        lista = listUFOV[y:y+5]
        maximo = max(lista)
        minimo = min(lista)
        dif2 = maximo-minimo
        #print(f"lista {lista[0]} {lista[1]} {lista[2]} {lista[3]} {lista[4]}\n")
       
        if dif2>dif: 
            dif = dif2
            lUFOVd = lista
            #print(f"max {lUFOVd[0]} {lUFOVd[1]} {lUFOVd[2]} {lUFOVd[3]} {lUFOVd[4]}\n")
       
        
 maximo = max(lUFOVd)
 minimo = min(lUFOVd)
 UFOVd = 100*(maximo - minimo)/(maximo + minimo)
 print(f"\n UFOV diferencial {UFOVd}")

 # CFOV diferencial

 dif = 0
 
 for y in range(0,len(listCFOV),5):
        lista2 = listCFOV[y:y+5]
        maximo = max(lista2)
        minimo = min(lista2)
        dif2 = maximo-minimo
        #print(f"lista {lista2[0]} {lista2[1]} {lista2[2]} {lista2[3]} {lista2[4]}\n")

        if dif2>dif: 
            dif = dif2
            lCFOVd = lista2
            #print(f"max {lCFOVd[0]} {lCFOVd[1]} {lCFOVd[2]} {lCFOVd[3]} {lCFOVd[4]}")
            

 maximo = max(lCFOVd)
 minimo = min(lCFOVd)
 CFOVd = 100*(maximo - minimo)/(maximo + minimo)
 print(f"\n CFOV diferencial {CFOVd}\n")

 return UFOVi, CFOVi, UFOVd, CFOVd


# Filas y columnas que queramos restar
# Cambiar los valores de int a float 
l=list(d1.shape)
filas = l[0]
columnas = l[1]
#df = np.empty([filas,columnas], dtype=float)
#df[:] = d1
print(f"\nEl archivo DICOM a analizar tiene:\n {filas} filas y \n {columnas} columnas.\n")
restaf = input("Introduzca el número de filas de píxels que desea restar:\n")
restac = input("Introduzca el número de columnas de píxels que desea restar:\n")
rf = int(restaf)
rc = int(restac)

matriz = suavizado(d1, rf//2, rc//2)
# print("filtro de suavizado")
# print(matriz)
# print("")

# anulamos pixeles con un número de cuentas anormalmente bajo
# def anular(matriz):
#         l = list(matriz.shape)
#         filas = l[0]
#         columnas = l[1]
#         restof=round(0.25*filas)
#         restoc=round(0.25*columnas)
        
#         rf = range(0, restof//2)
#         lf = list(rf)
#         rc = range(0, restoc//2)
#         lc = list(rc)
#         rf2 = range(filas-restof//2, filas)
#         lf2 = list(rf2)
#         rc2 = range(columnas-restoc//2, columnas)
#         lc2 = list(rc2)
        
#         mCFOV = np.delete(matriz, lf2, 0)
#         mCFOV = np.delete(mCFOV, lc2, 1)
#         mCFOV = np.delete(mCFOV, lf, 0)
#         mCFOV = np.delete(mCFOV, lc, 1)
        
#         # Restamos los valores nulos y calculamos la media de cuentas del CFOV
#         datos = mCFOV.flatten()
#         lista= list()
#         for x in range (0, len(datos)):
#          if datos[x] != 0:
#              lista.append(float(datos[x]))
    
#         media = np.mean(lista)
#         normaliz = matriz.flatten()
#         for z in range(0, len(normaliz)):
#             if float(normaliz[z]) < media:
#                normaliz[z] = 0

#         resultado = np.reshape(normaliz, (filas,columnas))
#         return resultado
# matriz = anular(d1)
# print("anulamos pixeles con un número de cuentas anormalmente bajo")
# print(matriz)
# print("")

# Anulamos vecinos de los valores nulos
def zeros(matriz):
 datos = np.copy(matriz)
 lm = list(matriz.shape)
 filasm = lm[0]
 columnasm = lm[1]
 n=0
 for i in range(0,filasm):
  for j in range(0, columnasm):
        if matriz[i,j] == 0:
          vecinosf = [i-1, i, i+1]
          vecinosc = [j-1 , j, j+1]
          for ni in vecinosf:
             for nj in vecinosc:
                 if 0 <= ni < filasm and 0 <= nj < columnasm:
                        datos[i][nj] = 0
                        datos[ni][j] = 0
 
 return datos
matriz = zeros(matriz)
print(matriz)
# print("anulamos vecinos de 0s")
# print(d)
# print("")


#Calculamos la uniformidad
uniformidad(matriz)
 


    
