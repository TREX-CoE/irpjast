#!/usr/bin/python

import numpy as np
import sys, getopt

from jastlib import getCoefList, get_sphere_distribution, generateElsAroundPoints

def main(argv):
   Natom = 2
   Ratio = 5
   inputfile = ''
   outputfile = ''
   try:
      opts, args = getopt.getopt(argv,"ha:r:",["atom=","ratio="])
   except getopt.GetoptError:
      print('test.py -a <natoms> -r <ratio>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print('test.py -a <inputfile> -r <outputfile>')
         sys.exit()
      elif opt in ("-a", "--atom"):
         Natom = int(arg)
      elif opt in ("-r", "--ratio"):
         Ratio = int(arg)


   Nelec = Ratio
   Nord = 5

   L1 = 20.0
   n = Natom # number of points
   dmin = 0.1 # min dist between points
   Ls = np.array([L1,L1,L1]) # lengths of the box
   shift = -10.0
   kappa = 2.0
   filename_atom = "geometry.txt"
   filename_coeffsa = "jast_coeffs.txt"
   (coeffsa, coeffsb, coeffsc) = getCoefList(Nord,n)
   coeffsall = np.concatenate((coeffsa,coeffsb,coeffsc))

   atomList,_,_ = get_sphere_distribution(n, dmin, Ls, maxiter=1e4, allow_wall=True)

   L1 = 15.0
   n = Nelec # number of points
   dmin = 0.1 # min dist between points
   Ls = np.array([L1,L1,L1]) # lengths of the box
   shift = -10.0
   kappa = 2.0
   filename_elec = "elec_coord.txt"

   rlist = generateElsAroundPoints(n,atomList,dmin)

   # Save file
   np.savetxt(filename_elec,rlist)
   np.savetxt(filename_atom,atomList)
   np.savetxt(filename_coeffsa,coeffsall)

if __name__ == "__main__":
   main(sys.argv[1:])
