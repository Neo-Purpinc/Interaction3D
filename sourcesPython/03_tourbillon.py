
import math
from Vector import *
import random

f = open("animations/03_tourbillon.pv","w")

time_step = 0.1
nb_points = 200
nb_frames = 600

max_dist = 3
# place 10 lignes et 10 colonnes de points
tous_les_points=[]
for point in range(nb_points):
    rayon = random.random()*max_dist
    # uniformally distributed in a sphere of radius y
    x = rayon * math.cos(2*math.pi/nb_points*point)
    z = rayon * math.sin(2*math.pi/nb_points*point)
    tous_les_points += [Vector(x,0,z)]

# ecriture de l'entete
f.write("#PV==\n" + str(time_step) + "\n" + str(nb_points) + "\n" + str(nb_frames) + "\n")
f.write("Point 0:199\n")
f.write("====\n")


#-------------- animation
# parcours des frames (images)
for fr in range(nb_frames):
    # parcours des points
    for n in range(nb_points):
        (x,y,z)=tous_les_points[n].get()
        rayon = math.sqrt(x*x + z*z)
        # make the points rotate around the center, with a speed proportional to their distance from the center
        x = rayon * math.cos(2*math.pi/nb_points*fr - rayon)
        z = rayon * math.sin(2*math.pi/nb_points*fr - rayon)
        f.write(str(x) + " " + str(y) + " " + str(z) + "\n")

f.close()
