#! /usr/bin/python

"""
Script to perform comparions between quadrature integrators 
and the analytic integrator. One can also do timing tests
here.
"""

import KCalc
from scipy import integrate
from scipy.special import spherical_jn as sphj
import time

#Functional form of the covariance integrand.
def BesselIntegrand(x,n,l,a,b):
 return x**n*sphj(l,a*x)*sphj(l,b*x)

#Testing values
n=6
l=4
a = 0.01
b = 50.0
x1 = 1e-3
x2 = 1e3

#Maximum number of quadrature subdivisions.
intlimit=1000

#Compute using both techniques. Compare results, times.
def Compare(x1,x2,l,n,a,b):
	start = time.clock()
	an1=KCalc.Calculate((x1,x2),l,n,a,b)
	end = time.clock()
	adelta = end-start
	start = time.clock()
	quad1= integrate.quad(BesselIntegrand,x1,x2,args=(n,l,a,b),limit=intlimit)
	end = time.clock()
	qdelta = end-start
	diff = abs(an1-quad1[0])
	quad_err = abs(quad1[1])
	print "Analytics: ","{:.10E}".format(float(an1))
	print "Quadrature:","{:.10E}".format(quad1[0])
	print "Diff: ", "{:.5E}".format(diff)
	print "QErr: ", "{:.5E}".format(quad_err)
	print
	print "Analytic Time: ",adelta," Seconds"
	print "Quadrature Time: ",qdelta," Seconds"
	print "Speedup: ",qdelta/adelta


#Loop	through values of l and execute comparisons.
for	l in range(4,10):
	print "n=",n
	print "l=",l
	Compare(x1,x2,l,n,a,b)
	print
	print

