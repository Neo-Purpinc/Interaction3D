
import math
from Vector import *
import random

f = open("animations/04_boom2.pv","w")

time_step = 0.1
nb_points = 600
nb_frames = 600
max_dist = 3
tous_les_points = []  

for i in range(nb_points):
    # place les points aléatoirement dans le plan y = 0
    x = random.uniform(-max_dist,max_dist)
    z = random.uniform(-max_dist,max_dist)
    tous_les_points.append(Vector(x,0,z))

# ecriture de l'entete
f.write("#PV==\n" + str(time_step) + "\n" + str(nb_points) + "\n" + str(nb_frames) + "\n")
f.write("Point 0:"+str(nb_points-1)+"\n")
f.write("====\n")

I = Vector(0,-1,0)

#-------------- animation
# parcours des frames (images)
for frame in range(nb_frames):
    for i in range(nb_points):
        (x,y,z)=tous_les_points[i].get()
        f.write(str(x) + " " + str(y) + " " + str(z) + "\n")
        vecteur_vitesse = tous_les_points[i] - I
        intensite = 1/vecteur_vitesse.module()
        x = x + intensite*vecteur_vitesse.x/100
        y = y + intensite*vecteur_vitesse.y/100
        z = z + intensite*vecteur_vitesse.z/100
        tous_les_points[i].set(x,y,z)

        


f.close()
