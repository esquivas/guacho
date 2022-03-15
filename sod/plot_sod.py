#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
from guacho_utils import *

#  ruta a donde estan los datos
path = '/home_diable/esquivel/guacho-dev/sod/BIN/'
nout = 20

#  lectura de los arreglos 3D scale=False es por que no usamos escalamientos
#  todo el calculo fue en unidades de código
rho = readbin3d_all(nout=nout,neq=0,path=path,verbose=False, scale=False )
vx  = readbin3d_all(nout=nout,neq=1,path=path,verbose=False, scale=False )
vy  = readbin3d_all(nout=nout,neq=2,path=path,verbose=False, scale=False )
vz  = readbin3d_all(nout=nout,neq=3,path=path,verbose=False, scale=False )
Pth = readbin3d_all(nout=nout,neq=4,path=path,verbose=False, scale=False )
Temp= Pth/rho

#  obtenemos los ejes (coordenadas de xi)
(x_axis, y_axis, z_axis) = get_axis(nout,path)


plt.ion()  #  modo interactivo, libera la terminal cuando aparece una figura

# Crear nueva figura de tamaño 10 x 8 "pulgadas" en una ventana con numero 1
plt.figure(figsize=(10,8), num = 1)

# Con subplot(filas, columnas, posición) especificamos que queremos un arreglo
# de sub-gráficas con ese número de filas y columnas, y que actualmente queremos
# graficar en esa posición (la numeración empieza en 1 y corre de izquierda a
# derecha, y de arriba hacia abajo)
plt.subplot(2,2,1)

# El comando plot() recibe un vector con los valores de x y otro con los valores
# de y como los primeros dos argumentos. Los demás argumentos son opciones,
# como el estilo de los marcadores y el color.
plt.plot(x_axis, rho[4,4,:], marker=".", color="blue")
# Título de esta (sub-)gráfica
plt.title("Densidad")
# Rango del eje x que queremos mostrar
plt.xlim(0, 1)
# Rango del eje y que queremos mostrar
plt.ylim(0, 1.05)


# Los siguientes 3 grupos grafican diferentes variables en subplots distintos
plt.subplot(2,2,2)
plt.plot(x_axis, vx[4,4,:], marker=".", color="green")
plt.title("Velocidad")
plt.xlim(0, 1)
plt.ylim(-0.05, 1.0)

plt.subplot(2,2,3)
plt.plot(x_axis, Pth[4,4,:], marker=".", color="darkorange")
plt.title("Presión")
plt.xlim(0, 1)
plt.ylim(0, 1.05)

plt.subplot(2,2,4)
plt.plot(x_axis, Temp[4,4,:], marker=".", color="red")
plt.title("Temperatura")
plt.xlim(0, 1)
plt.ylim(0.7, 1.2)

# Ajusta los márgenes entre subplots para que queden más estrechos pero sin
# sobreponer subplots
plt.tight_layout()

# Descomentar siguiente línea si se quiere guardar la figura a disco como imagen
# plt.savefig("figura.png")

# Muestra la figura en pantalla
plt.show()

