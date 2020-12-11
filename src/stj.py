#!/usr/bin/env python
# -*- coding: utf-8 -*-
from scipy.optimize import fsolve
from math import exp
from numpy import *




def solveStj(stj_type,T_fluid_in,A_ref_stj,theta,dir_irr,diff_irr,T_amb,stj_m_flow,stjFluidCp,delta_t_m_f):

        # --------------------------------------------------------------
        # Calculation of Solar Thermal Venetian Blinds (STJ) energy yield based on their performance curves [ISO9806]
        
        # Parameters: 
        # stj_type = Currently there are 3 STJ types, based on blinds extensions (0 for fully retracted and 1 for fully down), and slat angles| Accepted values ={1,2,3} {from config}
        # T_fluid  = Inlet fluid temperature [C] {from config}
        # A_ref = Collector's refrence area [m²] {from config}
        # Theta = Solar altitude angle (rad)
        # dir_irr = direct Irradiance beam [W/m²]
        # diff_irr = Diffuse Irradiance [W/m²]
        # T_amb = Ambient temperature [C]
        # stj_m_flow = collectors fluid masss flow [kg/s] {from config}
        # delta_t_m_f = mean temprature derivative  [C/s]
          
		if stj_type == 1:
			#CASE 1 Blinds Extensions = 1 (down) Beta = 82 deg (slat angles)
			def myFunction_1 (z):
				eta0 = 0.305
				alfa_1 = 3.05
				alfa_2 = 0.0
				alfa_5 = 0.138
				b_0 = 0.15
				k_b = 1 - b_0 * ((1/math.cos(theta*3.14/180))-1)
				k_d = 0.898
				Q_use_stj = z[0]
				T_m = z[1]
				F = empty((2))
				F[0] = (Q_use_stj/A_ref_stj) - ( eta0 * k_b * dir_irr + eta0 * k_d * diff_irr - alfa_1 * (T_m - T_amb) - alfa_2 * (T_m - T_amb)**2 - alfa_5 * delta_t_m_f) # 
				F[1] = Q_use_stj - ( stj_m_flow * stjFluidCp * ( 2*T_m - 2* T_fluid_in) )

				return F
			zGuess = array ([100,20])
			z = fsolve(myFunction_1,zGuess, xtol=1.49012e-08, factor=0.1)
			
		elif stj_type == 2:
			#CASE 2 Blinds Extensions = 0.5 (halfway down) Beta = 82 deg (slat angles)
			def myFunction_2 (z):
				eta0 = 0.222
				alfa_1 = 3.26
				alfa_2 = 0.0
				alfa_5 = 0.133
				b_0 = 0.16
				k_b = 1 - b_0 * ((1/math.cos(theta*3.14/180))-1)
				k_d = 0.898
				Q_use_stj = z[0]
				T_m = z[1]
				F = empty((2))
				F[0] = (Q_use_stj/A_ref_stj) - ( eta0 * k_b * dir_irr + eta0 * k_d * diff_irr - alfa_1 * (T_m - T_amb) - alfa_2 * (T_m - T_amb)**2 - alfa_5 * delta_t_m_f)
				F[1] = Q_use_stj - ( stj_m_flow * stjFluidCp * ( 2*T_m - 2* T_fluid_in) )

				return F
			zGuess = array ([100,20])
			z = fsolve(myFunction_2,zGuess, xtol=1.49012e-08, factor=0.1)
		elif stj_type == 3:
			#CASE 3 Blinds Extensions = 0.5 (halfway down) Beta = 45 deg (slat angles)
			def myFunction_3 (z):
				eta0 = 0.202
				alfa_1 = 2.67
				alfa_2 = 0.0
				alfa_5 = 0.201
				b_0 = -0.1
				k_b = 1 - b_0 * ((1/math.cos(theta*3.14/180))-1)
				k_d = 0.898 # needs to be fixed
				Q_use_stj = z[0]
				T_m = z[1]
				F = empty((2))
				F[0] = (Q_use_stj/A_ref_stj) - ( eta0 * k_b * dir_irr + eta0 * k_d * diff_irr - alfa_1 * (T_m - T_amb) - alfa_2 * (T_m - T_amb)**2 - alfa_5 * delta_t_m_f)
				F[1] = Q_use_stj - ( stj_m_flow * stjFluidCp * ( 2*T_m - 2* T_fluid_in) )
				return F
			zGuess = array ([100,20])
			z = fsolve(myFunction_3,zGuess, xtol=1.49012e-08, factor=0.1)
		heatflux = z[0]
		T_mean_f = z[1]
		heatFlux = max(heatflux,0)	
			
		return heatFlux,T_mean_f
