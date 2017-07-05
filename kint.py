#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Evaluates indefinite integrals of the form
\int x^n j_l(ax) j_l(bx) dx
for n=2, 4, 6 with a != b and 0 <= l <= 10 by hard-coded expressions
Warning: do not use at x = 0!
"""

from scipy import integrate
from mpmath import mp
import numpy as np
from math import cos, sin
from scipy.special import spherical_jn as sphj
import time

mp.dps=200

class Polynomial(object) :
		"""
		Computes polynomials of the form
		coeff[0] a^2n + coeff[1] a^2(n-1) b^2 + ...
		... + coeff[n-2] a^2 b^2(n-1) + coeff[n-1] b^2n
		a and b are set at initialization
		coefficients are passed in as a list
		"""

		def __init__(self, a, b) :
				self.a = a
				self.b = b
				self.a2 = a*a
				self.b2 = b*b
				self.stored = {}

		def get_list(self, n) :
				if n in self.stored :
						return self.stored[n]
				self.stored[n] = self._make_list(n)
				return self.stored[n]

		def _make_list(self, n) :
				"""
				Creates a list of monomials:
				a^2n, a^2(n-1) b^2, ..., a^2 b^2(n-1), b^2n
				"""
				# Creates a list of (a^2n, a^2(n-1), ..., a^2, 1)
				apoly = np.full(n, self.a2)
				apoly[0] = mp.mpf(1)
				apoly = np.flipud(np.cumprod(apoly))
				# Same for b, just reversed
				bpoly = np.full(n, self.b2)
				bpoly[0] = mp.mpf(1)
				bpoly = np.cumprod(bpoly)
				# Multiply the two together to get the list of monomials
				return apoly * bpoly

		def eval_poly(self, coeff) :
				"""
				Evaluates a polynomial of the form
				coeff[0] a^2n + coeff[1] a^2(n-1) b^2 + ...
				... + coeff[n-2] a^2 b^2(n-1) + coeff[n-1] b^2n
				"""
				return np.sum(coeff * self.get_list(len(coeff)))

class Term(object) :
		"""
		Describes a term of the form
		C x^n (ab)^m poly(a, b)
		where poly is a polynomial in a^2 and b^2
		"""
		def __init__(self, c, n, m, coeffs=[1]) :
			self.c = mp.mpf(c)
			self.n = mp.mpf(n)
			self.m = mp.mpf(m)
			self.coeffs = []
			for coeff in coeffs:
				self.coeffs.append(mp.mpf(coeff))

		def evaluate(self, x, poly) :
				"""Evaluates the term, given x and a poly object that stores a and b"""
				xn = x ** self.n
				abm = mp.mpf(1)
				if self.m != 0 :
						abm = (poly.a * poly.b) ** self.m
				return self.c * xn * abm * poly.eval_poly(self.coeffs)

class Part(object) :
		"""
		Describes a part of an integral as a sum of terms
		"""
		def __init__(self, terms) :
				"""Pass in a list of terms"""
				self.terms = terms

		def evaluate(self, x, poly) :
				"""Evaluates the part, given x and a poly object that stores a and b"""
				return sum([term.evaluate(x, poly) for term in self.terms])

class KnlInt(object) :
		"""
		Describes an integral as an appropriate sum of parts
		"""
		def __init__(self, n, l, csum, cdiff, ssum, sdiff) :
				"""Describe the integral, and provide a list of parts"""
				self.n = mp.mpf(n)
				self.l = mp.mpf(l)
				self.parts = [csum, cdiff, ssum, sdiff]

		def evaluate(self, x, poly) :
				"""Evaluates the integral, given x and a poly object that stores a and b"""
				# Evaluate each part (csum, cdiff, ssum, sdiff)
				csum, cdiff, ssum, sdiff = [part.evaluate(x, poly) for part in self.parts]
				# Compute the required coefficients
				abl = (poly.a * poly.b) ** (self.l + mp.mpf(1))
				apb = poly.a + poly.b
				amb = poly.a - poly.b
				# Evaluate each term
				cpterm = mp.cos(x * apb) * (csum + cdiff) / apb ** (self.n - mp.mpf(2))
				cmterm = mp.cos(x * amb) * (csum - cdiff) / amb ** (self.n - mp.mpf(2))
				spterm = mp.sin(x * apb) * (ssum + sdiff) / apb ** (self.n - mp.mpf(1))
				smterm = mp.sin(x * amb) * (ssum - sdiff) / amb ** (self.n - mp.mpf(1))
				# Compute result and check precision
				termlist = [cpterm,cmterm,spterm,smterm]
				termlist = [mp.mpf('0.25') / abl * term for term in termlist]
				# Compute the order of sum of abs(terms)
				abslist = [mp.fabs(term) for term in termlist]
				maxdigits = mp.floor(mp.log10(sum(abslist)))
				# Compute the order of abs(sum(terms))
				result = sum(termlist)
				resdigits = mp.floor(mp.log10(mp.fabs(result)))
				# Compute precision loss
				prec_loss = maxdigits-resdigits
				prec_remaining = mp.dps - prec_loss
				print "Precision remaining",prec_remaining
				# Demand that at least machine precision remains FIXME: more?
				if prec_remaining<13:
					sys.exit("Insufficient precision")
				# Add the lot together
				# print mp.mpf('0.25') / abl * (cpterm + cmterm + spterm + smterm)
				return result

def K2Int(l, x, poly) :
		"""
		Computes the integral
		\int x^2 j_l(ax) j_l(bx) dx
		which is known analytically
		"""
		a = poly.a
		b = poly.b
		coeff = x*x/(a*a-b*b)
		if l == 0 :
				terms = sphj(l, a*x)*cos(b*x)/x - cos(a*x)/x*sphj(l, b*x)
		else :
				terms = b*sphj(l, a*x)*sphj(l-1, b*x) - a*sphj(l-1, a*x)*sphj(l, b*x)
		return coeff*terms

integrals = {}

# n=4, l=0
n = 4
l = 0
csum = Part([])
cdiff = Part([Term(-4, 1, 0)])
ssum = Part([Term(-4, 2, 1)])
sdiff = Part([Term(4, 0, 0), Term(-2, 2, 0, [1, 1])])
integrals[(n, l)] = KnlInt(n, l, csum, cdiff, ssum, sdiff)

# n=4, l=1
n = 4
l = 1
csum = Part([Term(8, 1, 1)])
cdiff = Part([Term(2, 1, 0, [1, 1])])
ssum = Part([Term(-12, 0, 1), Term(2, 2, 1, [1, 1])])
sdiff = Part([Term(4, 2, 2), Term(-4, 0, 0, [1, 1])])
integrals[(n, l)] = KnlInt(n, l, csum, cdiff, ssum, sdiff)

# n=4, l=2
n = 4
l = 2
csum = Part([Term(36, -1, 1), Term(-6, 1, 1, [1, 1])])
cdiff = Part([Term(-16, 1, 2), Term(18, -1, 0, [1, 1])])
ssum = Part([Term(-4, 2, 3), Term(36, 0, 1, [1, 1])])
sdiff = Part([Term(-2, 2, 2, [1, 1]), Term(2, 0, 0, [3, 32, 3])])
integrals[(n, l)] = KnlInt(n, l, csum, cdiff, ssum, sdiff)

# n=4, l=3
n = 4
l = 3
csum = Part([Term(300, -3, 1), Term(28, 1, 3), Term(-210, -1, 1, [1, 1])])
cdiff = Part([Term(150, -3, 0, [1, 1]), Term(12, 1, 2, [1, 1]), Term(-30, -1, 0, [1, 12, 1])])
ssum = Part([Term(600, -2, 1, [1, 1]), Term(2, 2, 3, [1, 1]), Term(-2, 0, 1, [15, 116, 15])])
sdiff = Part([Term(4, 2, 4), Term(-144, 0, 2, [1, 1]), Term(150, -2, 0, [1, 6, 1])])
integrals[(n, l)] = KnlInt(n, l, csum, cdiff, ssum, sdiff)

# n=4, l=4
n = 4
l = 4
csum = Part([Term(8820, -5, 1), Term(-7770, -3, 1, [1,1]), Term(-20, 1, 3, [1,1]), Term(30, -1, 1, [7, 60, 7])])
cdiff = Part([Term(-44, 1, 4), Term(4410, -5, 0, [1, 1]), Term(1110, -1, 2, [1,1]), Term(-420, -3, 0, [4, 29, 4])])
ssum = Part([Term(-4, 2, 5), Term(17640, -4, 1, [1, 1]), Term(400, 0, 3, [1,1]), Term(-210, -2, 1, [11, 50, 11])])
sdiff = Part([Term(-2, 2, 4, [1, 1]), Term(4410, -4, 0, [1, 6, 1]), Term(-210, -2, 0, [1, 35, 35, 1]), Term(6, 0, 2, [15, 104, 15])])
integrals[(n, l)] = KnlInt(n, l, csum, cdiff, ssum, sdiff)

# n=4, l=5
n = 4
l = 5
csum = Part([Term(510300, -7, 1), Term(64, 1, 5), Term(-470610, -5, 1, [1,1]), Term(-3990, -1, 3, [1,1]), Term(420, -3, 1, [63, 326, 63])])
cdiff = Part([Term(255150, -7, 0, [1,1]), Term(30, 1, 4, [1,1]), Term(-420, -1, 2, [2, 15, 2]), Term(210, -3, 0, [9, 443, 443, 9]), Term(-5670, -5, 0, [19, 128, 19])])
ssum = Part([Term(1020600, -6, 1, [1,1]), Term(2, 2, 5, [1,1]), Term(210, -2, 1, [9, 215, 215, 9]), Term(-5670, -4, 1, [31, 122, 31]), Term(-2, 0, 3, [105, 692, 105])])
sdiff = Part([Term(4, 2, 6), Term(-900, 0, 4, [1,1]), Term(255150, -6, 0, [1,6,1]), Term(-22680, -4, 0, [1, 22, 22, 1]), Term(420, -2, 2, [37, 150, 37])])
integrals[(n, l)] = KnlInt(n, l, csum, cdiff, ssum, sdiff)

# n=4, l=6
n = 4
l = 6
csum = Part([Term(48024900, -9, 1), Term(-45218250, -7, 1, [1, 1]), Term(-42, 1, 5, [1, 1]), Term(2520, -1, 3, [1, 7, 1]), Term(-630, -3, 1, [33, 995, 995, 33]), Term(1890, -5, 1, [1749, 7712, 1749])])
cdiff = Part([Term(-88, 1, 6), Term(24012450, -9, 0, [1, 1]), Term(11340, -1, 4, [1, 1]), Term(-623700, -7, 0, [17, 111, 17]), Term(-1260, -3, 2, [159, 710, 159]), Term(1890, -5, 0, [187, 5418, 5418, 187])])
ssum = Part([Term(-4, 2, 7), Term(96049800, -8, 1, [1, 1]), Term(1764, 0, 5, [1, 1]), Term(-6300, -2, 3, [11, 42, 11]), Term(7560, -4, 1, [55, 754, 754, 55]), Term(-103950, -6, 1, [177, 662, 177])])
sdiff = Part([Term(-2, 2, 6, [1, 1]), Term(24012450, -8, 0, [1, 6, 1]), Term(-3150, -2, 2, [3, 61, 61, 3]), Term(-103950, -6, 0, [25, 483, 483, 25]), Term(4, 0, 4, [105, 673, 105]), Term(1890, -4, 0, [11, 1205, 4040, 1205, 11])])
integrals[(n, l)] = KnlInt(n, l, csum, cdiff, ssum, sdiff)

# n=4, l=7
n = 4
l = 7
csum = Part([Term(6640533900, -11, 1), Term(116, 1, 7), Term(-6328372050, -9, 1, [1, 1]), Term(-27468, -1, 5, [1, 1]), Term(3118500, -7, 1, [169, 692, 169]), Term(-20790, -5, 1, [325, 5446, 5446, 325]), Term(2520, -3, 3, [407, 1675, 407])])
cdiff = Part([Term(3320266950, -11, 0, [1, 1]), Term(56, 1, 6, [1, 1]), Term(-252, -1, 4, [25, 168, 25]), Term(-28378350, -9, 0, [53, 340, 53]), Term(1260, -3, 2, [99, 2390, 2390, 99]), Term(623700, -7, 0, [104, 2471, 2471, 104]), Term(-20790, -5, 0, [13, 2032, 7452, 2032, 13])])
ssum = Part([Term(13281067800, -10, 1, [1, 1]), Term(2, 2, 7, [1, 1]), Term(3150, -2, 3, [11, 205, 205, 11]), Term(-12, 0, 5, [63, 397, 63]), Term(-28378350, -8, 1, [95, 346, 95]), Term(20790, -6, 1, [4017, 44359, 44359, 4017]), Term(-20790, -4, 1, [13, 831, 2560, 831, 13])])
sdiff = Part([Term(4, 2, 8), Term(-3136, 0, 6, [1, 1]), Term(3320266950, -10, 0, [1, 6, 1]), Term(-56756700, -8, 0, [7, 127, 127, 7]), Term(12600, -2, 4, [19, 70, 19]), Term(-83160, -4, 2, [44, 487, 487, 44]), Term(20790, -6, 0, [299, 18932, 58290, 18932, 299])])
integrals[(n, l)] = KnlInt(n, l, csum, cdiff, ssum, sdiff)

# n=4, l=8
n = 4
l = 8
csum = Part([Term(1264255492500, -13, 1), Term(-1214451488250, -11, 1, [1, 1]), Term(-72, 1, 7, [1, 1]), Term(1260, -1, 5, [11, 72, 11]), Term(-13860, -3, 3, [39, 830, 830, 39]), Term(-40540500, -7, 1, [48, 625, 625, 48]), Term(28378350, -9, 1, [3855, 15124, 3855]), Term(20790, -5, 1, [195, 16224, 53692, 16224, 195])])
cdiff = Part([Term(-148, 1, 8), Term(632127746250, -13, 0, [1, 1]), Term(59220, -1, 6, [1, 1]), Term(-7662154500, -11, 0, [38, 241, 38]), Term(-27720, -3, 4, [147, 575, 147]), Term(28378350, -9, 0, [510, 10907, 10907, 510]), Term(20790, -5, 2, [3107, 40158, 40158, 3107]), Term(-8108100, -7, 0, [15, 1259, 4182, 1259, 15])])
ssum = Part([Term(-4, 2, 9), Term(2528510985000, -12, 1, [1, 1]), Term(5184, 0, 7, [1, 1]), Term(-138600, -2, 5, [5, 18, 5]), Term(-1351350, -6, 1, [99, 3396, 9410, 3396, 99]), Term(-3831077250, -10, 1, [139, 498, 139]), Term(83160, -4, 3, [260, 2549, 2549, 260]), Term(56756700, -8, 1, [345, 3409, 3409, 345])])
sdiff = Part([Term(-2, 2, 8, [1, 1]), Term(632127746250, -12, 0, [1, 6, 1]), Term(-34650, -2, 4, [3, 53, 53, 3]), Term(-3831077250, -10, 0, [21, 367, 367, 21]), Term(4, 0, 6, [315, 1963, 315]), Term(-1351350, -6, 0, [3, 862, 7335, 7335, 862, 3]), Term(28378350, -8, 0, [60, 3017, 8862, 3017, 60]), Term(20790, -4, 2, [91, 4525, 13240, 4525, 91])])
integrals[(n, l)] = KnlInt(n, l, csum, cdiff, ssum, sdiff)

# n=4, l=9
n = 4
l = 9
csum = Part([Term(316653859021500, -15, 1), Term(184, 1, 9), Term(-305907687335250, -13, 1, [1, 1]), Term(-116820, -1, 7, [1, 1]), Term(13860, -3, 5, [975, 3686, 975]), Term(-270270, -5, 3, [1555, 17262, 17262, 1555]), Term(2554051500, -11, 1, [11373, 43390, 11373]), Term(-28378350, -9, 1, [21998, 249171, 249171, 21998]), Term(40540500, -7, 1, [68, 2903, 8565, 2903, 68])])
cdiff = Part([Term(158326929510750, -15, 0, [1, 1]), Term(90, 1, 8, [1, 1]), Term(-3960, -1, 6, [7, 45, 7]), Term(-716411445750, -13, 0, [103, 648, 103]), Term(6930, -3, 4, [273, 5363, 5363, 273]), Term(1277025750, -11, 0, [3145, 62991, 62991, 3145]), Term(4054050, -7, 0, [17, 6765, 65753, 65753, 6765, 17]), Term(-270270, -5, 2, [120, 7329, 22736, 7329, 120]), Term(-56756700, -9, 0, [833, 52514, 164475, 52514, 833])])
ssum = Part([Term(633307718043000, -14, 1, [1, 1]), Term(2, 2, 9, [1, 1]), Term(6930, -2, 5, [39, 665, 665, 39]), Term(-4, 0, 7, [495, 3061, 495]), Term(-238803815250, -12, 1, [573, 2030, 573]), Term(255405150, -10, 1, [21947, 202709, 202709, 21947]), Term(-1891890, -4, 3, [5, 215, 608, 215, 5]), Term(1351350, -6, 1, [51, 7586, 55503, 55503, 7586, 51]), Term(-4054050, -8, 1, [13600, 350695, 913962, 350695, 13600])])
sdiff = Part([Term(4, 2, 10), Term(-8100, 0, 8, [1, 1]), Term(158326929510750, -14, 0, [1, 6, 1]), Term(-7567560, -4, 4, [13, 118, 118, 13]), Term(-955215261000, -12, 0, [22, 375, 375, 22]), Term(13860, -2, 6, [127, 450, 127]), Term(-8108100, -8, 0, [323, 48975, 361340, 361340, 48975, 323]), Term(510810300, -10, 0, [1037, 45881, 130820, 45881, 1037]), Term(1351350, -6, 2, [1062, 27015, 70126, 27015, 1062])])
integrals[(n, l)] = KnlInt(n, l, csum, cdiff, ssum, sdiff)

# n=4, l=10
n = 4
l = 10
csum = Part([Term(100863567447142500, -17, 1), Term(-97855355786438250, -15, 1, [1, 1]), Term(-110, 1, 9, [1, 1]), Term(1980, -1, 7, [26, 165, 26]), Term(-12870, -3, 5, [441, 8203, 8203, 441]), Term(-21709437750, -11, 1, [10811, 112109, 112109, 10811]), Term(716411445750, -13, 1, [13471, 50392, 13471]), Term(-4054050, -7, 1, [323, 62135, 507299, 507299, 62135, 323]), Term(270270, -5, 3, [680, 34517, 102480, 34517, 680]), Term(16216200, -9, 1, [93347, 2888861, 7948878, 2888861, 93347])])
cdiff = Part([Term(-224, 1, 10), Term(50431783723571250, -17, 0, [1, 1]), Term(214830, -1, 8, [1, 1]), Term(-353907254200500, -15, 0, [67, 419, 67]), Term(-25740, -3, 6, [1519, 5606, 1519]), Term(716411445750, -13, 0, [1919, 36748, 36748, 1919]), Term(270270, -5, 4, [7820, 78617, 78617, 7820]), Term(-43418875500, -11, 0, [456, 24339, 73330, 24339, 456]), Term(-8108100, -7, 2, [3910, 118660, 324617, 118660, 3910]), Term(4054050, -9, 0, [15181, 3001979, 24809428, 24809428, 3001979, 15181])])
ssum = Part([Term(-4, 2, 11), Term(201727134894285000, -16, 1, [1, 1]), Term(12100, 0, 9, [1, 1]), Term(7567560, -4, 5, [49, 422, 422, 49]), Term(-25740, -2, 7, [157, 550, 157]), Term(-176953627100250, -14, 1, [251, 882, 251]), Term(955215261000, -12, 1, [2052, 18113, 18113, 2052]), Term(-1351350, -6, 3, [7718, 164231, 408030, 164231, 7718]), Term(8108100, -8, 1, [8075, 602055, 3737318, 3737318, 602055, 8075]), Term(-1240539300, -10, 1, [19323, 423039, 1059772, 423039, 19323])])
sdiff = Part([Term(-2, 2, 10, [1, 1]), Term(50431783723571250, -16, 0, [1, 6, 1]), Term(-176953627100250, -14, 0, [39, 653, 653, 39]), Term(-12870, -2, 6, [49, 815, 815, 49]), Term(6, 0, 8, [495, 3044, 495]), Term(1891890, -4, 4, [20, 783, 2162, 783, 20]), Term(-1351350, -6, 2, [459, 49179, 326326, 326326, 49179, 459]), Term(238803815250, -12, 0, [817, 33313, 93060, 33313, 817]), Term(-620269650, -10, 0, [2242, 250543, 1691711, 1691711, 250543, 2242]), Term(4054050, -8, 0, [323, 208318, 3842665, 9287180, 3842665, 208318, 323])])
integrals[(n, l)] = KnlInt(n, l, csum, cdiff, ssum, sdiff)

# n=6, l=0
n = 6
l = 0
csum = Part([Term(-16, 3, 1)])
cdiff = Part([Term(48, 1, 0), Term(-8, 3, 0, [1, 1])])
ssum = Part([Term(48, 2, 1), Term(-8, 4, 1, [1, 1])])
sdiff = Part([Term(-48, 0, 0), Term(24, 2, 0, [1, 1]), Term(-2, 4, 0, [1, 6, 1])])
integrals[(n, l)] = KnlInt(n, l, csum, cdiff, ssum, sdiff)

# n=6, l=1
n = 6
l = 1
csum = Part([Term(-80, 1, 1), Term(16, 3, 1, [1, 1])])
cdiff = Part([Term(-16, 1, 0, [1, 1]), Term(2, 3, 0, [1, 14, 1])])
ssum = Part([Term(80, 0, 1), Term(-56, 2, 1, [1, 1]), Term(2, 4, 1, [1, 6, 1])])
sdiff = Part([Term(16, 0, 0, [1, 1]), Term(8, 4, 2, [1, 1]), Term(-8, 2, 0, [1, 12, 1])])
integrals[(n, l)] = KnlInt(n, l, csum, cdiff, ssum, sdiff)

# n=6, l=2
n = 6
l = 2
csum = Part([Term(168, 1, 1, [1, 1]), Term(-2, 3, 1, [3, 26, 3])])
cdiff = Part([Term(-32, 3, 2, [1, 1]), Term(6, 1, 0, [5, 54, 5])])
ssum = Part([Term(-240, 0, 1, [1, 1]), Term(-8, 4, 3, [1, 1]), Term(12, 2, 1, [5, 26, 5])])
sdiff = Part([Term(-2, 4, 2, [1, 6, 1]), Term(-48, 0, 0, [1, 9, 1]), Term(6, 2, 0, [1, 35, 35, 1])])
integrals[(n, l)] = KnlInt(n, l, csum, cdiff, ssum, sdiff)

# n=6, l=3
n = 6
l = 3
csum = Part([Term(1800, -1, 1, [1, 1]), Term(56, 3, 3, [1, 1]), Term(-30, 1, 1, [11, 58, 11])])
cdiff = Part([Term(450, -1, 0, [1, 6, 1]), Term(4, 3, 2, [3, 22, 3]), Term(-6, 1, 0, [5, 191, 191, 5])])
ssum = Part([Term(2, 4, 3, [1, 6, 1]), Term(-6, 2, 1, [5, 111, 111, 5]), Term(60, 0, 1, [25, 98, 25])])
sdiff = Part([Term(8, 4, 4, [1, 1]), Term(-12, 2, 2, [19, 78, 19]), Term(6, 0, 0, [35, 701, 701, 35])])
integrals[(n, l)] = KnlInt(n, l, csum, cdiff, ssum, sdiff)

# n=6, l=4
n = 6
l = 4
csum = Part([Term(29400, -3, 1, [1, 1]), Term(70, 1, 1, [3, 73, 73, 3]), Term(-4, 3, 3, [5, 34, 5]), Term(-1050, -1, 1, [15, 58, 15])])
cdiff = Part([Term(-88, 3, 4, [1, 1]), Term(7350, -3, 0, [1, 6, 1]), Term(-2100, -1, 0, [1, 21, 21, 1]), Term(2, 1, 2, [855, 3634, 855])])
ssum = Part([Term(-8, 4, 5, [1, 1]), Term(14700, -2, 1, [3, 10, 3]), Term(-50, 0, 1, [63, 737, 737, 63]), Term(4, 2, 3, [155, 582, 155])])
sdiff = Part([Term(-2, 4, 4, [1, 6, 1]), Term(7350, -2, 0, [1, 15, 15, 1]), Term(2, 2, 2, [45, 847, 847, 45]), Term(-2, 0, 0, [105, 7710, 24394, 7710, 105])])
integrals[(n, l)] = KnlInt(n, l, csum, cdiff, ssum, sdiff)

# n=6, l=5
n = 6
l = 5
csum = Part([Term(1428840, -5, 1, [1, 1]), Term(128, 3, 5, [1, 1]), Term(1260, -1, 1, [27, 349, 349, 27]), Term(-13230, -3, 1, [71, 250, 71]), Term(-6, 1, 3, [1015, 3938, 1015])])
cdiff = Part([Term(357210, -5, 0, [1, 6, 1]), Term(-120, 1, 2, [7, 142, 142, 7]), Term(-13230, -3, 0, [11, 185, 185, 11]), Term(2, 3, 4, [15, 98, 15]), Term(630, -1, 0, [3, 284, 930, 284, 3])])
ssum = Part([Term(714420, -4, 1, [3, 10, 3]), Term(2, 4, 5, [1, 6, 1]), Term(-42, 2, 3, [5, 87, 87, 5]), Term(-13230, -2, 1, [21, 187, 187, 21]), Term(6, 0, 1, [315, 14980, 43418, 14980, 315])])
sdiff = Part([Term(8, 4, 6, [1, 1]), Term(357210, -4, 0, [1, 15, 15, 1]), Term(-276, 2, 4, [5, 18, 5]), Term(3000, 0, 2, [7, 67, 67, 7]), Term(-26460, -2, 0, [1, 43, 120, 43, 1])])
integrals[(n, l)] = KnlInt(n, l, csum, cdiff, ssum, sdiff)

# n=6, l=6
n = 6
l = 6
csum = Part([Term(123492600, -7, 1, [1, 1]), Term(168, 1, 3, [15, 277, 277, 15]), Term(-2, 3, 5, [21, 134, 21]), Term(-561330, -5, 1, [151, 522, 151]), Term(1890, -3, 1, [2519, 25109, 25109, 2519]), Term(-1890, -1, 1, [11, 628, 1890, 628, 11])])
cdiff = Part([Term(-176, 3, 6, [1, 1]), Term(30873150, -7, 0, [1, 6, 1]), Term(-2245320, -5, 0, [6, 97, 97, 6]), Term(-11340, -1, 2, [23, 241, 241, 23]), Term(60, 1, 4, [287, 1062, 287]), Term(1890, -3, 0, [209, 11109, 32620, 11109, 209])])
ssum = Part([Term(-8, 4, 7, [1, 1]), Term(61746300, -6, 1, [3, 10, 3]), Term(48, 2, 5, [56, 197, 56]), Term(-187110, -4, 1, [157, 1267, 1267, 157]), Term(-168, 0, 3, [555, 4807, 4807, 555]), Term(3780, -2, 1, [132, 3341, 8526, 3341, 132])])
sdiff = Part([Term(-2, 4, 6, [1, 6, 1]), Term(30873150, -6, 0, [1, 15, 15, 1]), Term(12, 2, 4, [35, 583, 583, 35]), Term(1890, -2, 0, [11, 1902, 13559, 13559, 1902, 11]), Term(-187110, -4, 0, [17, 602, 1610, 602, 17]), Term(-30, 0, 2, [315, 12460, 34506, 12460, 315])])
integrals[(n, l)] = KnlInt(n, l, csum, cdiff, ssum, sdiff)

# n=6, l=7
n = 6
l = 7
csum = Part([Term(16232416200, -9, 1, [1, 1]), Term(232, 3, 7, [1, 1]), Term(-44594550, -7, 1, [255, 874, 255]), Term(1260, -1, 3, [1067, 9969, 9969, 1067]), Term(166320, -5, 1, [4667, 41660, 41660, 4667]), Term(-4, 1, 5, [10395, 37406, 10395]), Term(-6930, -3, 1, [1131, 33919, 90468, 33919, 1131])])
cdiff = Part([Term(4058104050, -9, 0, [1, 6, 1]), Term(8, 3, 6, [7, 44, 7]), Term(-44594550, -7, 0, [41, 651, 651, 41]), Term(-28, 1, 4, [225, 3931, 3931, 225]), Term(-6930, -3, 0, [39, 9069, 71176, 71176, 9069, 39]), Term(1260, -1, 2, [99, 4497, 12880, 4497, 99]), Term(41580, -5, 0, [1807, 76435, 214132, 76435, 1807])])
ssum = Part([Term(8116208100, -8, 1, [3, 10, 3]), Term(2, 4, 7, [1, 6, 1]), Term(-4, 2, 5, [189, 3065, 3065, 189]), Term(-14864850, -6, 1, [283, 2197, 2197, 283]), Term(-6930, -2, 1, [39, 3946, 24927, 24927, 3946, 39]), Term(20790, -4, 1, [5239, 101496, 243810, 101496, 5239]), Term(2, 0, 3, [17325, 620424, 1672366, 620424, 17325])])
sdiff = Part([Term(8, 4, 8, [1, 1]), Term(4058104050, -8, 0, [1, 15, 15, 1]), Term(-8, 2, 6, [595, 2064, 595]), Term(140, 0, 4, [2295, 18761, 18761, 2295]), Term(-29729700, -6, 0, [16, 529, 1390, 529, 16]), Term(-13860, -2, 2, [321, 6393, 15484, 6393, 321]), Term(20790, -4, 0, [325, 31671, 196644, 196644, 31671, 325])])
integrals[(n, l)] = KnlInt(n, l, csum, cdiff, ssum, sdiff)

# n=6, l=8
n = 6
l = 8
csum = Part([Term(2988240255000, -11, 1, [1, 1]), Term(-8, 3, 7, [9, 56, 9]), Term(-5533778250, -9, 1, [383, 1306, 383]), Term(36, 1, 5, [385, 6491, 6491, 385]), Term(12162150, -7, 1, [13215, 111469, 111469, 13215]), Term(-41580, -1, 3, [13, 519, 1440, 519, 13]), Term(62370, -3, 1, [65, 8203, 56168, 56168, 8203, 65]), Term(-1621620, -5, 1, [1495, 33659, 84308, 33659, 1495])])
cdiff = Part([Term(-296, 3, 8, [1, 1]), Term(747060063750, -11, 0, [1, 6, 1]), Term(-11067556500, -9, 0, [31, 487, 487, 31]), Term(-41580, -1, 4, [129, 1123, 1123, 129]), Term(60, 1, 6, [1491, 5270, 1491]), Term(-6486480, -5, 0, [20, 2479, 16828, 16828, 2479, 20]), Term(12162150, -7, 0, [1370, 52121, 142386, 52121, 1370]), Term(20790, -3, 2, [3679, 83999, 211260, 83999, 3679])])
ssum = Part([Term(-8, 4, 9, [1, 1]), Term(1494120127500, -10, 1, [3, 10, 3]), Term(-5533778250, -8, 1, [147, 1117, 1117, 147]), Term(24, 2, 7, [327, 1124, 327]), Term(-36, 0, 5, [25795, 203261, 203261, 25795]), Term(291060, -2, 3, [91, 1583, 3708, 1583, 91]), Term(-810810, -4, 1, [185, 9707, 52092, 52092, 9707, 185]), Term(8108100, -6, 1, [3230, 54919, 127750, 54919, 3230])])
sdiff = Part([Term(-2, 4, 8, [1, 6, 1]), Term(747060063750, -10, 0, [1, 15, 15, 1]), Term(84, 2, 6, [15, 239, 239, 15]), Term(145530, -2, 2, [13, 1026, 6017, 6017, 1026, 13]), Term(-5533778250, -8, 0, [17, 542, 1410, 542, 17]), Term(4054050, -6, 0, [470, 35923, 207655, 207655, 35923, 470]), Term(-30, 0, 4, [3465, 116760, 309286, 116760, 3465]), Term(-810810, -4, 0, [5, 1937, 28227, 63630, 28227, 1937, 5])])
integrals[(n, l)] = KnlInt(n, l, csum, cdiff, ssum, sdiff)

# n=6, l=9
n = 6
l = 9
csum = Part([Term(730739674665000, -13, 1, [1, 1]), Term(368, 3, 9, [1, 1]), Term(41580, -1, 5, [429, 3563, 3563, 429]), Term(-976924698750, -11, 1, [535, 1818, 535]), Term(-12, 1, 7, [14685, 51274, 14685]), Term(851350500, -9, 1, [49759, 405409, 405409, 49759]), Term(3243240, -5, 1, [935, 59405, 344101, 344101, 59405, 935]), Term(-270270, -3, 3, [1865, 36289, 87780, 36289, 1865]), Term(-60810750, -7, 1, [13158, 254303, 613806, 254303, 13158])])
cdiff = Part([Term(182684918666250, -13, 0, [1, 6, 1]), Term(-2520, 1, 6, [11, 181, 181, 11]), Term(2, 3, 8, [45, 278, 45]), Term(-976924698750, -11, 0, [87, 1357, 1357, 87]), Term(-810810, -3, 2, [40, 3753, 23555, 23555, 3753, 40]), Term(20790, -1, 4, [91, 3348, 9090, 3348, 91]), Term(-121621500, -7, 0, [425, 39561, 247196, 247196, 39561, 425]), Term(425675250, -9, 0, [10727, 383712, 1031794, 383712, 10727]), Term(810810, -5, 0, [85, 43930, 725636, 1696226, 725636, 43930, 85])])
ssum = Part([Term(365369837332500, -12, 1, [3, 10, 3]), Term(2, 4, 9, [1, 6, 1]), Term(-12, 2, 7, [165, 2597, 2597, 165]), Term(-325641566250, -10, 1, [631, 4729, 4729, 631]), Term(-1891890, -2, 3, [5, 342, 1909, 1909, 342, 5]), Term(-20270250, -6, 1, [3230, 125051, 611383, 611383, 125051, 3230]), Term(425675250, -8, 1, [17697, 279340, 637158, 279340, 17697]), Term(6, 0, 5, [45045, 1457940, 3816038, 1457940, 45045]), Term(810810, -4, 1, [85, 17395, 214701, 464310, 214701, 17395, 85])])
sdiff = Part([Term(8, 4, 10, [1, 1]), Term(182684918666250, -12, 0, [1, 15, 15, 1]), Term(-48, 2, 8, [255, 871, 255]), Term(1800, 0, 6, [1309, 10061, 10061, 1309]), Term(-3783780, -2, 4, [32, 511, 1170, 511, 32]), Term(-651283132500, -10, 0, [37, 1153, 2980, 1153, 37]), Term(851350500, -8, 0, [697, 46965, 260146, 260146, 46965, 697]), Term(810810, -4, 2, [2020, 78421, 383895, 383895, 78421, 2020]), Term(-40540500, -6, 0, [68, 13886, 171043, 369670, 171043, 13886, 68])])
integrals[(n, l)] = KnlInt(n, l, csum, cdiff, ssum, sdiff)

# n=6, l=10
n = 6
l = 10
csum = Part([Term(228624086213523000, -15, 1, [1, 1]), Term(-2, 3, 9, [55, 338, 55]), Term(440, 1, 7, [117, 1892, 1892, 117]), Term(-231400896977250, -13, 1, [711, 2410, 711]), Term(151966064250, -11, 1, [91903, 731341, 731341, 91903]), Term(-12870, -1, 5, [441, 15348, 41030, 15348, 441]), Term(90090, -3, 3, [2040, 160491, 949025, 949025, 160491, 2040]), Term(8108100, -7, 1, [216087, 9860255, 51638802, 51638802, 9860255, 216087]), Term(-482431950, -9, 1, [630515, 11102120, 26156106, 11102120, 630515]), Term(-270270, -5, 1, [4845, 1243210, 17052484, 38054674, 17052484, 1243210, 4845])])
cdiff = Part([Term(-448, 3, 10, [1, 1]), Term(57156021553380750, -15, 0, [1, 6, 1]), Term(-925603587909000, -13, 0, [29, 450, 450, 29]), Term(-25740, -1, 6, [2009, 16143, 16143, 2009]), Term(2, 1, 8, [161865, 560254, 161865]), Term(-964863900, -9, 0, [22534, 1787775, 10595035, 10595035, 1787775, 22534]), Term(90090, -3, 4, [28380, 497591, 1171170, 497591, 28380]), Term(-1081080, -5, 2, [32895, 1493435, 7805639, 7805639, 1493435, 32895]), Term(21709437750, -11, 0, [71117, 2442851, 6497480, 2442851, 71117]), Term(4054050, -7, 0, [15827, 4095912, 56365545, 125906008, 56365545, 4095912, 15827])])
ssum = Part([Term(-8, 4, 11, [1, 1]), Term(114312043106761500, -14, 1, [3, 10, 3]), Term(-77133632325750, -12, 1, [853, 6331, 6331, 853]), Term(-2200, 0, 7, [2457, 18553, 18553, 2457]), Term(4, 2, 9, [4565, 15522, 4565]), Term(1261260, -2, 5, [363, 5469, 12320, 5469, 363]), Term(-270270, -4, 3, [44540, 1445229, 6646695, 6646695, 1445229, 44540]), Term(-137837700, -8, 1, [211983, 6922955, 31917318, 31917318, 6922955, 211983]), Term(8683775100, -10, 1, [302860, 4556895, 10257402, 4556895, 302860]), Term(2702700, -6, 1, [26163, 2653768, 26982649, 55545304, 26982649, 2653768, 26163])])
sdiff = Part([Term(-2, 4, 10, [1, 6, 1]), Term(57156021553380750, -14, 0, [1, 15, 15, 1]), Term(2, 2, 8, [1485, 23167, 23167, 1485]), Term(630630, -2, 4, [60, 3743, 20181, 20181, 3743, 60]), Term(-77133632325750, -12, 0, [101, 3098, 7970, 3098, 101]), Term(4341887550, -10, 0, [50027, 3121140, 16805745, 16805745, 3121140, 50027]), Term(-2, 0, 6, [315315, 9921780, 25747834, 9921780, 315315]), Term(1351350, -6, 0, [969, 769029, 21161479, 92938987, 92938987, 21161479, 769029, 969]), Term(-270270, -4, 2, [2295, 342890, 3799599, 7983360, 3799599, 342890, 2295]), Term(-68918850, -8, 0, [21698, 3276813, 36459210, 76693582, 36459210, 3276813, 21698])])
integrals[(n, l)] = KnlInt(n, l, csum, cdiff, ssum, sdiff)

# Testing


def BesselIntegrand(x,n,l,a,b):
	return x**n*sphj(l,a*x)*sphj(l,b*x)


n=6
l=4
a = mp.mpf('0.01')
b = mp.mpf('50.0')
x1 = mp.mpf('1e-3')
x2 = mp.mpf('1e3')

def Compare(x1,x2,n,l,a,b):
	start = time.clock()
	poly = Polynomial(a, b)
	result1 = integrals[(n, l)].evaluate(x1, poly)
	result2 = integrals[(n, l)].evaluate(x2, poly)
	an1 = result2-result1
	end = time.clock()
	adelta = end-start


	start = time.clock()
	quad1= integrate.quad(BesselIntegrand,float(x1),float(x2),args=(n,l,float(a),float(b)),limit=1000)
	end = time.clock()
	qdelta = end-start
	diff = abs(float(an1)-quad1[0])
	quad_err = abs(quad1[1])
	print "Analytics: ","{:.10E}".format(float(an1))
	print "Quadrature:","{:.10E}".format(quad1[0])
	print "Diff: ", "{:.5E}".format(diff)
	print "QErr: ", "{:.5E}".format(quad_err)
	print
	print "Analytic Time: ",adelta," Seconds"
	print "Quadrature Time: ",qdelta," Seconds"
	print "Speedup: ",qdelta/adelta

for l in range(4,10):
	print "n=",n
	print "l=",l
	Compare(x1,x2,n,l,a,b)
	print
	print

#print(K2Int(l, x, mypoly))
