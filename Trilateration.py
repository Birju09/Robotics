#\usr\bin\env python

## This program uses Least squares method to find out the co-ordinates x,y,z from the given measurements for Localization
## Since there are three unknowns and more than 3 equations this is a redundant system of equations.
## The equations are obtained as :
## Let P(x,y,z) be the point whose distance from the satellites (beacon) is given
## Hence 'u' number of sphere equations will be formed which need to be solved, in the given data there are four sphere equations
## So Let first sphere be : (x-a)^2 + (y-a)^2 + (z-a)^2 = r1^2
## Similary the second equation will be (x-u)^2 + (y-v)^2 + (z-w)^2 = r2^2
## Subtracting these two equations will give us the equation -(a-u)x + -(b-v)y + -(c-w)z = (r1^2 - r2^2) - (a^2-u^2) - (b^2-v^2) - (c^2-w^2)
## Similarly four equations can be formed for the data given by subtracting pairs of the sphere equations
## Then these four equations are solved for three unknowns using the least square formulas for Ax = b as x = (A^TA)^-1 A^T b
## for matrix calculations numpy has been used
## The answers are within the accuracy of 10e-5
## The program also works for more then 4 measurements, we have checked for upto 24 measurements, but it can work for more measurments also

import sys
import math
from numpy import *


def trilaterate3D(distances):
    xl = []
    yl = []
    zl = []
    dl = []
    
    for a in range(len(distances)):
        xl.append(distances[a][0])
        yl.append(distances[a][1])
        zl.append(distances[a][2])
        dl.append(distances[a][3])
    [u,v] = [len(distances),len(distances[0])]
    A = empty((u,3),float)
    b1  = []
    ran = (math.factorial(u)/(math.factorial(u-2)*math.factorial(2)))  ##calculating the number of row of the A matrix, i.e. using nCr to find number of two pairs which can be subtracted to create the equations
    d = 1
    while d <= (ran-2):
        b =0
        while (b+d)<u:
            row = array([[-xl[b%u] + xl[(b+d)%u],-yl[b%u] + yl[(b+d)%u],-zl[b%u] + zl[(b+d)%u]]])                  ##row creation for A
            br = [((pow(dl[b%u],2)-pow(dl[(b+d)%u],2)) - (pow(xl[b%u],2)-pow(xl[(b+d)%u],2)) - (pow(yl[b%u],2)-pow(yl[(b+d)%u],2)) - (pow(zl[b%u],2)-pow(zl[(b+d)%u],2)))/2]     ##row creation for b
            A = append(A,row,axis=0)
            b1.append(br)
            b = b + 1
        d = d + 1
        
    A = delete(A,[range(u)],0)
    b = array(b1)
    Att = dot((A.T),A)
    X =  dot(dot(linalg.inv(Att),A.T),b)                     ##X = (A^TA)^-1A^Tb
    
    return [X[0][0],X[1][0],X[2][0],1] 

if __name__ == "__main__":

    # Retrive file name for input data
    if (len(sys.argv) == 1):
        print "Please enter data file name."
        exit()

    filename = sys.argv[1]

    # Read data
    lines = [line.rstrip('\n') for line in open(filename)]
    distances = []
    for line in range(0, len(lines)):
        distances.append(map(float, lines[line].split(' ')))

    # Print out the data
    print "The input four points and distances, in the format of [x, y, z, d], are:"
    for p in range(0, len(distances)):
        print distances[p]

        # Call the function and compute the location
    location = trilaterate3D(distances)
    print
    print "The location of the point is: " + str(location)
