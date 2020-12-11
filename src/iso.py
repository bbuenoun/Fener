#!/usr/bin/env python
# -*- coding: utf-8 -*-


from numpy import ones, zeros
import ISO15099_v31 as iso_modules


def compute(
    win,
    surf,
    frame,
    tempOutAir,
    tempInAir,
    extCHTC,
    CHTCwin,
    skyTemp,
    vfWinSurf,
    vfWinWin,
    vfWinFrame,
    irrWinLay,
    dynamicFlag=True,
):

    for win_id, one_win in enumerate(win):

        # --------------------------------------------------------------
        # Calculation of the windows energy transfers according to ISO15099:

        # First, we initialize the parameters and inputs:

        # Calculation of E_gv_int, the total longwave radiative heat flux coming from the room's surface towards the window

        E_gv_int = 0
        if dynamicFlag:
            # Important hypothesis: the emissivity of the surface doesn't appear here (factor 1),
            # because we add the radiation coming from the room and reflected from the elements towards the window
            # and assume this part is near to (1-Epsilon)*iso_modules.SB*(T_element**4)
            for surf_id, one_surf in enumerate(surf):
                E_gv_int += (
                    vfWinSurf[win_id * len(surf) + surf_id]
                    * iso_modules.SB
                    * (one_surf.temp[0] ** 4)
                )  # Part coming from the surfaces
            for win_id2, one_win2 in enumerate(win):
                E_gv_int += (
                    vfWinWin[win_id * len(win) + win_id2]
                    * iso_modules.SB
                    * (one_win2.intTemp ** 4)
                )  # Part coming from the other windows
            for frame_id, one_frame in enumerate(frame):
                E_gv_int += (
                    vfWinFrame[win_id * len(frame) + frame_id]
                    * iso_modules.SB
                    * (one_frame.intTemp ** 4)
                )  # Part coming from the frames

        # XXX eventually read the inputs from weather data for the HTCs:

        IRtrans = zeros(win[win_id].numPane)  # Hypothesis: no IR transmittance
        frontRefl = (
            ones(win[win_id].numPane) - win[win_id].frontEmis
        )  # Hypothesis: no IR transmittance
        backRefl = (
            ones(win[win_id].numPane) - win[win_id].backEmis
        )  # Hypothesis: no IR transmittance
        # airflow_rates=XXX
        # print(IRtrans)
        parameters_and_inputs = iso_modules.Inputs(
            # PARAMETERS AND INPUTS:
            win[win_id].numPane,  # Number of solid layers
            [-77777, -77777] + win[win_id].airFrac[1:].tolist(),  # air_ratio
            [-77777, -77777] + win[win_id].argFrac[1:].tolist(),  # argon_ratio
            [-77777, -77777, 0, 0, 0],  # krypton_ratio XXX can be implemented
            [-77777, -77777, 0, 0, 0],  # xenon_ratio XXX can be implemented
            [-77777]
            + win[win_id].paneCond.tolist()
            + [-77777],  # [W/(m*K)] lambda : thermal conductivity
            [-77777]
            + win[win_id].paneThick.tolist()
            + [-77777],  # [m] Thickness of the solid layers
            [-77777, -77777]
            + win[win_id].gapThick[1:].tolist(),  # [m] Thickness of the air layers
            [-77777]
            + win[win_id].frontEmis.tolist()
            + [-77777],  # [-] Longwave emissivity of the front of the solid layer
            [-77777]
            + win[win_id].backEmis.tolist()
            + [-77777],  # [-] Longwave emissivity of the back of the solid layers
            [-77777]
            + IRtrans.tolist()
            + [-77777],  # [-] Longwave transmissivity of solid layers
            [-77777]
            + frontRefl.tolist()
            + [-77777],  # [-] Longwave reflectivity of the front of the solid layers
            [-77777]
            + backRefl.tolist()
            + [-77777],  # [-] Longwave reflectivity of the back of the solid layers
            win[win_id].height,  # [m] the height of the glazed area
            win[win_id].length,  # [m] the width of the glazed area
            # INPUTS:
            [-7777]
            + irrWinLay[win_id, 0 : win[win_id].numPane].tolist()
            + [
                -7777
            ],  # XXX attention when number of layer varies      # [W/m²] Energy absorbed by each layer
            [
                -77777,
                -77777,
                0,
                0,
                0,
                -77777,
            ],  # [m³/h] Airflow rate of the corresponding  air gaps. airflow_rates[i] correspond to the airflow betweeen layers i and i-1.
            tempOutAir,  # [K] exterior temperature
            tempInAir,  # [K] room temperature
            [
                -77777,
                -77777,
                0,
                0,
                0,
                -77777,
            ],  # [K] Temperature of the inlet air in the gap
            extCHTC,  # [W/(m²*K)] Exterior convective heat transfer coeficient
            CHTCwin[win_id],  # [W/(m²*K)] Interior convective heat transfer coeficient
            0.5 * iso_modules.SB * tempOutAir ** 4
            + 0.5
            * skyTemp,  # XXXiso_modules.SB*T_rm_ex**4,    # [W/m²] Longwave irradiance incoming at external glass portion from the exterior
            E_gv_int,  # [W/m²] Longwave irradiance incoming at internal glass portion from the interior
        )
        the_outputs = iso_modules.Outputs()
        if dynamicFlag:
            # print(irrWinLay[win_id,0:win[win_id].numPane,it].tolist())
            # Then, we run the simulation and store the outputs:
            iso_modules.calculate_ISO_15099(parameters_and_inputs, the_outputs)
            # if irrWinLay[win_id,0:win[win_id].numPane,it].all()>0: iso_modules.print_outputs(parameters_and_inputs,the_outputs)
            # -------------------------------------------------------------

            # transmission of window outer surface temperature:
            win[win_id].extTemp = the_outputs.T_ft[1]
            # transmission of window inner surface temperature:
            win[win_id].intTemp = the_outputs.T_b[-2]
            # only for stj
            win[win_id].blindTemp = the_outputs.T_ft[2]
            win[win_id].qInt = the_outputs.q_int
        else:
            # Then, we run the simulation and store the outputs:
            iso_modules.calculate_ISO_15099(
                parameters_and_inputs, the_outputs, 0, CHTCwin[win_id], extCHTC
            )
            # print("the_outputs.q_int"+str(the_outputs.q_int))
            # secondary heat transmission to interior
            win[win_id].qInt = the_outputs.q_int
            # print("win[win_id].qInt IN ISO.PY: "+str(win[win_id].qInt))
            # U-value
            win[win_id].uValue = the_outputs.U_gv

    return win
