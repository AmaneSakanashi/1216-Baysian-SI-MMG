import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
# import ray

from shipsim.step import ShipManeuver
from my_src.ddcma import *
from my_src.read_bound import *
from my_src.read_train import *
from .obj_function import Obj_function 


ship = ShipManeuver(dt_sim=0.1, solve_method="euler")

def deg2rad(deg):
    return deg * np.pi / 180

class SI:
    def __init__(self):
        self.dt_sim = 0.1
        
        self.no_files, self.no_timestep,\
        self.set_action_train, self.set_state_train, self.set_wind_train = Read_train.read_csv(self)

        self.Dim = 57
        self.N = 2*self.Dim
        
        bound = Set_bound()   
        self.LOWER_BOUND, self.UPPER_BOUND, self.FLAG_PERIODIC, self.period_length = bound.set_param_bound(self.N)

        ### initial state ###
        ### w_xxx : Weight of the Obj. term
        self.w_max = 1e+2
        self.w_noise = [0.03,0.01,0.03,0.01, deg2rad(0.1), deg2rad(0.1)]
        self.w_pen = 1e+6
        self.w_overflow = [0.3,0.1,0.3,0.1, deg2rad(1), deg2rad(1)]
        self.const_Obj = self.N * (math.log( 2 * math.pi + 1 ))

        self.FLAG_OVERFLOW = False

    # @ray.remote
    def calc_step(self, x):
            obj = Obj_function(Dim=57)
            w_noise, w_max, w_overflow, FLAG_OVERFLOW  = self.w_noise, self.w_max, self.w_overflow, self.FLAG_OVERFLOW
            update_params = x.copy()        

            # t_list = []  
            for j in range(self.no_files):
                func = 0

                state_train = self.set_state_train[j,:,:]
                action_train = self.set_action_train[j,:,:]
                wind_train = self.set_wind_train[j,:,:]
                # ---------------------------------------------------------------------------------------------------------------
                ### calculating trj. loop ###
                for i in range(int(self.no_timestep)):
                # Reset to training data value every 100s
                    if (i%(100/self.dt_sim)==0):
                        init_state = state_train[i,:]
                        self.state_sim = init_state
                        self.t = 0.0

                    t_n, state_sim_n = ship.step(self.t, update_params, self.state_sim,  action_train[i], wind_train[i])
            #obj_func
                    func_i, FLAG_OVERFLOW = obj.J_norm(state_sim_n, state_train[i], w_noise, w_max, w_overflow, FLAG_OVERFLOW )
                    func += (-1)*func_i
            # update
                    self.t, self.state_sim = t_n, state_sim_n

                if FLAG_OVERFLOW == True:
                     func = func * self.w_pen

                func_lkh = obj.J_lkh(update_params)
                func += (-1)* func_lkh

            return func

    def fobj(self, x):
            params_list = mirror(x, self.LOWER_BOUND, self.UPPER_BOUND, self.FLAG_PERIODIC)
            cand_results = []
            
            for set_params in params_list:
                
                ## -------------------------------------------------------

                m_params = set_params[0:self.Dim]
                v_params = set_params[self.Dim:]

                ## For degug ---------------------------------------------
                # mean_result = pd.read_csv("log/Opt_mean_result.csv", header=None)
                # var_result = pd.read_csv("log/Opt_var_result.csv", header=None)
                # m_params = np.array(mean_result.values.flatten())
                # v_params = np.array(var_result.values.flatten())
                ## -------------------------------------------------------
                det_v_params = abs(np.sum(v_params))

                results = []
                gen_params = np.array([np.random.normal(m, v, 2) 
                        for m, v in zip(m_params,v_params)]).T
                for i in range(gen_params.ndim):
                    # per_result_id = actor.trj_culc.remote(actor,gen_params[i]) 
                    per_result_id = self.calc_step(gen_params[i]) 
                    results.append(per_result_id)

                # results = ray.get(results)
                results_all = np.mean(results)   + \
                                0.5 *( math.log(det_v_params)+ self.const_Obj)
                cand_results.append(-1 * results_all)

            return  cand_results



