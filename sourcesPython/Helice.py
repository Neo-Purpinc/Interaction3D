import math
from Vector import *

#-------------------------------------------------------
# Classe representant des helices definies par :
# x = rayon * cos(pulsation * t)
# y = pente * t
# z = rayon * sin(pulsation * t)
# ou t est le parametre de l'helice. Pour t=0, on obtient
# un point a la base de l'helice (en y=0) et, au fur et
# a mesure que t augmente, on obtient des points plus
# haut (si pente>0) sur l'helice et pour t negatif,
# on obtient des points plus bas sur l'helice.

class Helice:

    #-------
    # Constructeur

    def __init__(self, rayon, pulsation, pente):
        self.rayon = rayon
        self.pulsation = pulsation
        self.pente = pente

    #-------
    # t doit etre un flottant et represente le parametre
    # Valeur de retour : un Vecteur representant la position
    # du point de parametre t sur l'helice. En parcourant
    # toutes les valeurs de t, on obtient tous les points
    # de l'helice.

    def courbe(self, t):
        x = self.rayon * math.cos(self.pulsation*t)
        y = self.pente * t
        z = self.rayon * math.sin(self.pulsation*t)
        return Vector(x,y,z)

    #------
    # t doit etre un flottant et represente le parametre
    # Valeur de retour : un Vecteur representant le vecteur
    # tangent a la courbe au point de parametre t sur l'helice.

    def tangente(self, t):
        x = - self.rayon * self.pulsation * math.sin(self.pulsation*t)
        y = self.pente
        z = self.rayon * self.pulsation * math.cos(self.pulsation*t)
        return Vector(x,y,z)

    #------
    # y doit etre un flottant representant une hauteur dans l'espace
    # Valeur de retour : un vecteur representant la position du point
    # de l'helice qui se trouve a la hauteur y.

    def pt_a_la_hauteur(self, y):
        return self.courbe(y/self.pente)

    #------
    # y doit etre un flottant representant une hauteur dans l'espace
    # Valeur de retour : un vecteur representant la tangente a l'helice
    # au point de hauteur y.

    def tang_a_la_hauteur(self, y):
        return self.tangente(y/self.pente)

    #------
    # p doit etre un vecteur representant la position d'un point p
    # Valeur de retour : un flottant representant la distance
    # euclidienne entre le point p et le point de l'helice aui se
    # trouve a la meme hauteur que p.

    def distance_avec_point(self, p):
        ph = self.pt_a_la_hauteur(p.y)
        return p.distanceAvec(ph)
