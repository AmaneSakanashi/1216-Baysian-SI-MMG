import numpy as np
import pandas as pd
import mmg_fort as mf

class OptParam:
    def __init__(self):
        self.mean_result = pd.read_csv("opt_result/Opt_mean_result.csv", header=None)
        self.var_result = pd.read_csv("opt_result/Opt_var_result.csv", header=None)

        self.output_path = "log/sim_data/"
        principal_particulars = pd.read_csv(
            'inputfiles/principal_particulars_EssoOsaka3m.csv', header=0, index_col=0
        )
        parameter_init = pd.read_csv(
            "inputfiles/MMG_params_EssoOsaka3m.csv", header=0, index_col=0
        )
        #
        self.pp_vector = np.empty(len(principal_particulars))
        ### Principal Particulars ###
        self.pp_vector[0] = principal_particulars.at["lpp", "value"]
        self.pp_vector[1] = principal_particulars.at["breadth", "value"]
        self.pp_vector[2] = principal_particulars.at["draft", "value"]
        self.pp_vector[3] = principal_particulars.at["mass_nd", "value"]
        self.pp_vector[4] = principal_particulars.at["x_lcg", "value"]
        # Propeller
        self.pp_vector[5] = principal_particulars.at["dia_prop", "value"]
        self.pp_vector[6] = principal_particulars.at["pitchr", "value"]
        # Rudder
        self.pp_vector[7] = principal_particulars.at["area_rudder", "value"]
        self.pp_vector[8] = principal_particulars.at["lambda_rudder", "value"]
        self.pp_vector[9] = principal_particulars.at["x_location_rudder", "value"]
        # parameters for wind force computation
        self.pp_vector[10] = principal_particulars.at[
            "area_projected_trans", "value"
        ]  # AT
        self.pp_vector[11] = principal_particulars.at[
            "area_projected_lateral", "value"
        ]  # AL
        self.pp_vector[12] = principal_particulars.at["AOD", "value"]  # AOD
        self.pp_vector[13] = principal_particulars.at[
            "LCW", "value"
        ]  # C or LCW, midship to center of AL
        self.pp_vector[14] = principal_particulars.at[
            "LCBR", "value"
        ]  # CBR, midship to center of AOD(superstructure or bridge)
        self.pp_vector[15] = principal_particulars.at[
            "HBR", "value"
        ]  #  Height from free surface to top of the superstructure (bridge) (m)
        self.pp_vector[16] = principal_particulars.at[
            "HC", "value"
        ]  # Hc, hight of center of lateral projected area
        self.pp_vector[17] = principal_particulars.at[
            "SBW", "value"
        ]  # ship breadth fro wind force computation
        self.pp_vector[18] = principal_particulars.at[
            "Lz", "value"
        ]  # Lz --- Acting position of sway force from center of gravity (it is necessary to calculate roll moment)
        
        #
        self.mmg_params_vector = np.empty(len(parameter_init))
        ### Parameter Init ###
        self.mmg_params_vector[0] = parameter_init.at["massx_nd", "value"]
        self.mmg_params_vector[1] = parameter_init.at["massy_nd", "value"]
        self.mmg_params_vector[2] = parameter_init.at["IzzJzz_nd", "value"]
        # Hull
        self.mmg_params_vector[3] = parameter_init.at["xuu_nd", "value"]
        self.mmg_params_vector[4] = parameter_init.at["xvr_nd", "value"]
        self.mmg_params_vector[5] = parameter_init.at["yv_nd", "value"]
        self.mmg_params_vector[6] = parameter_init.at["yr_nd", "value"]
        self.mmg_params_vector[7] = parameter_init.at["nv_nd", "value"]
        self.mmg_params_vector[8] = parameter_init.at["nr_nd", "value"]
        self.mmg_params_vector[9] = parameter_init.at["coeff_drag_sway", "value"]
        self.mmg_params_vector[10] = parameter_init.at["cry_cross_flow", "value"]
        self.mmg_params_vector[11] = parameter_init.at["crn_cross_flow", "value"]
        self.mmg_params_vector[12] = parameter_init.at["coeff_drag_aft", "value"]
        # Propeller
        self.mmg_params_vector[13] = parameter_init.at["t_prop", "value"]
        self.mmg_params_vector[14] = parameter_init.at["w_prop_zero", "value"]
        self.mmg_params_vector[15] = parameter_init.at["tau_prop", "value"]
        self.mmg_params_vector[16] = parameter_init.at["coeff_cp_prop", "value"]
        self.mmg_params_vector[17] = parameter_init.at["xp_prop_nd", "value"]
        self.mmg_params_vector[18:21] = np.array(
            [
                parameter_init.at["kt_coeff0", "value"],
                parameter_init.at["kt_coeff1", "value"],
                parameter_init.at["kt_coeff2", "value"],
            ]
        )
        self.mmg_params_vector[21:29] = np.array(
            [
                parameter_init.at["ai_coeff_prop0", "value"],
                parameter_init.at["ai_coeff_prop1", "value"],
                parameter_init.at["ai_coeff_prop2", "value"],
                parameter_init.at["ai_coeff_prop3", "value"],
                parameter_init.at["ai_coeff_prop4", "value"],
                parameter_init.at["ai_coeff_prop5", "value"],
                parameter_init.at["ai_coeff_prop6", "value"],
                parameter_init.at["ai_coeff_prop7", "value"],
            ]
        )
        self.mmg_params_vector[29:37] = np.array(
            [
                parameter_init.at["bi_coeff_prop0", "value"],
                parameter_init.at["bi_coeff_prop1", "value"],
                parameter_init.at["bi_coeff_prop2", "value"],
                parameter_init.at["bi_coeff_prop3", "value"],
                parameter_init.at["bi_coeff_prop4", "value"],
                parameter_init.at["bi_coeff_prop5", "value"],
                parameter_init.at["bi_coeff_prop6", "value"],
                parameter_init.at["bi_coeff_prop7", "value"],
            ]
        )
        self.mmg_params_vector[37:41] = np.array(
            [
                parameter_init.at["ci_coeff_prop0", "value"],
                parameter_init.at["ci_coeff_prop1", "value"],
                parameter_init.at["ci_coeff_prop2", "value"],
                parameter_init.at["ci_coeff_prop3", "value"],
            ]
        )

        # Rudder
        self.mmg_params_vector[41] = parameter_init.at["t_rudder", "value"]
        self.mmg_params_vector[42] = parameter_init.at["ah_rudder", "value"]
        self.mmg_params_vector[43] = parameter_init.at["xh_rudder_nd", "value"]
        self.mmg_params_vector[44] = parameter_init.at["kx_rudder", "value"]
        self.mmg_params_vector[45] = parameter_init.at["epsilon_rudder", "value"]
        self.mmg_params_vector[46] = parameter_init.at["lr_rudder_nd", "value"]
        self.mmg_params_vector[47] = parameter_init.at["gammaN_rudder", "value"]
        self.mmg_params_vector[48] = parameter_init.at["gammaP_rudder", "value"]
        self.mmg_params_vector[49] = parameter_init.at["kx_rudder_reverse", "value"]
        self.mmg_params_vector[50] = parameter_init.at["cpr_rudder", "value"]
        self.mmg_params_vector[51] = parameter_init.at["XX0", "value"]
        self.mmg_params_vector[52] = parameter_init.at["XX1", "value"]
        self.mmg_params_vector[53] = parameter_init.at["XX3", "value"]
        self.mmg_params_vector[54] = parameter_init.at["XX5", "value"]
        self.mmg_params_vector[55] = parameter_init.at["YY1", "value"]
        self.mmg_params_vector[56] = parameter_init.at["YY3", "value"]
        self.mmg_params_vector[57] = parameter_init.at["YY5", "value"]
        self.mmg_params_vector[58] = parameter_init.at["NN1", "value"]
        self.mmg_params_vector[59] = parameter_init.at["NN2", "value"]
        self.mmg_params_vector[60] = parameter_init.at["NN3", "value"]

    def update(self):           
            m_params = np.array(self.mean_result.values.flatten())
            v_params = np.array(self.var_result.values.flatten())  

            update_params = np.random.normal(m_params, v_params)

            update_index = np.array([
                        0,1,2, # principal particular
                        4,5,6,7,8,9,10,11,12, # Hull force parameters
                        13,14,15,16,17, # Propeller
                        21,22,23,24,25,26,27,28, # Coeffs for Yp
                        29,30,31,32,33,34,35,36, # Coeffs for Np
                        37,38,39,40, # Coeffs for Xp
                        41,42,43,44,45,46,47,48,49,50, # Rudder
                        51,52,53,54,55,56,57,58,59,60 # Wind force coefficients
                                 ]) 
            update_mmg_params = self.mmg_params_vector

            for i in range(len(update_index)):
                update_mmg_params[update_index[i]] = update_params[i]

            return update_mmg_params
    
    def ode_rhs(self,x, delta_rudder, n_prop, Wind_velocity, Wind_Direction, update_params):
        # [delta_rudder, n_prop] = u
        # [Wind_velocity, Wind_Direction] = w
        n_bt = 0
        n_st = 0
        # print("MMG_pamams : ", self.mmg_params_vector)
        dx = mf.esso_osaka_vector_input.mmg_lowspeed_model(
            x,
            delta_rudder,
            n_prop,
            n_bt,
            n_st,
            Wind_Direction,
            Wind_velocity,
            self.pp_vector,
            update_params,
        )
        return dx