
import math
from Vector import *

f = open("animations/02_cercle.pv","w")

time_step = 0.5
nb_points = 1
nb_frames = 800

# ecriture de l'entete
f.write("#PV==\n" + str(time_step) + "\n" + str(nb_points) + "\n" + str(nb_frames) + "\n")
f.write("Point 0:0\n")
f.write("====\n")


#-------------- animation
rayon = 1
vitesse_rotation = 0.1
position_depart = Vector(0,0,0)

# parcours des frames (images)
for frame in range(nb_frames):
    (x,y,z)=position_depart.get()
    x=rayon*math.cos(vitesse_rotation*frame)
    z=rayon*math.sin(vitesse_rotation*frame)
    f.write(str(x) + " " + str(y) + " " + str(z) + "\n")

f.close()
