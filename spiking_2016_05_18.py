# -*- coding: utf-8 -*-

""" spiking_2016_05_18.py

    Model of Cortical Spiking Neurons as proposed by Izhikevich (2003)

    Copyright 2016 Tri-Peter Shrive

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Tri-Peter Shrive
    Anklamer Str 13
    10115 Berlin
    Deutschland

    +49 17 62087558
    Tri.Shrive@gmail.com

"""

import numpy as np
import matplotlib.pyplot as plt


def rk4(a, b, c, d, I, v0, u0, t0, Tt, n):

	Fv = [0] * (n + 1)
	Fu = [0] * (n + 1)
	Ft = [0] * (n + 1)

	h  = (Tt - t0) / float(n)

	Fv[0] = v = v0
	Fu[0] = u = u0


	for i in range(1, n + 1):
		k1_Dv, k1_Du = Dv(v, u, I, c), Du(v, u, a, b, d),

		k2_Dv, k2_Du = Dv(v + h * 0.5 * k1_Dv, u + h * 0.5 * k1_Du, I, c), Du(v + h * 0.5 * k1_Dv, u + h * 0.5 * k1_Du, a, b, d)

		k3_Dv, k3_Du = Dv(v + h * 0.5 * k2_Dv, u + h * 0.5 * k2_Du, I, c), Du(v + h * 0.5 * k2_Dv, u + h * 0.5 * k2_Du, a, b, d)

		k4_Dv, k4_Du = Dv(v + h * k3_Dv, u + h * k3_Du, I, c), Du(v + h * k3_Dv, u + h * k3_Du, a, b, d)

		Ft[i] = t = t0 + i * h
		Fv[i] = v = v + h * (k1_Dv + k2_Dv + k2_Dv + k3_Dv + k3_Dv + k4_Dv) / 6
		if Fv[i] >= 30:
			Fv[i] = 30
		Fu[i] = u = u + h * (k1_Du + k2_Du + k2_Du + k3_Du + k3_Du + k4_Du) / 6

		if v >= 30:
			v = c
			u = u + d

	return Fv, Fu, Ft

def Dv(v, u, I, c):
	return 0.04 * v**2 + 5 * v + 140 - u + I

def Du(v, u, a, b, d):
	return a * (b * v - u)

def overtime(Fv, Fu, Ft, t0, Tt):

	save = np.array([ Fv, Fu, Ft ])
	np.savetxt("spiking.csv", save, delimiter = ",")

	n = 20000

	Fv, Fu, Ft = Fv[0:n], Fu[0:n], Ft[0:n]

	fig = plt.figure(1)
	fig.canvas.set_window_title('')
	fig.suptitle('Spiking Model by Izhikevich (2003)', fontsize=12, fontweight='bold')

	av = fig.add_subplot(211)
	av.set_ylabel('v(t)')
	plt.plot( Ft, Fv, 'r' )

	au = fig.add_subplot(212)
	au.set_ylabel('u(t)')
	au.set_xlabel('t')
	plt.plot( Ft, Fu, 'r' )
	plt.show()

	return

if __name__ == '__main__':

	v0, u0 = -56, -110
	a, b, c, d, I = 0.2, 2, -56, -16, -99

	t0 , Tt = 0, 3000
	n =  20 * Tt # RK4 step size
	
	Fv, Fu, Ft = rk4(a, b, c, d, I, v0, u0, t0, Tt, n)

	overtime(Fv, Fu, Ft, t0, Tt)

