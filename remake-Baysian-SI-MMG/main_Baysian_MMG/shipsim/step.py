import numpy as np
from .ode_rhs import *

MMG = EssoOsakaMMG()

class ShipManeuver:
    def __init__(self,dt_sim, solve_method):
            self.solve_method = solve_method
            self.dt_sim = dt_sim

    def step(self, t, update_params, ship_state, ship_action, ship_wind):
        """simulation step

        Args:
            
            last_logging (bool, optional): If you want to save last (next) state, change this to True. Defaults to False.

        Returns:
            ObsType: Observation variable for the next step
            bool: Handle to determine the end of simulation
            dict: Additional infomation
        """
        # t, ship_state = self.t, self.ship_state

        # start loop
        ### time ###
        t_n = t + self.dt_sim
        ### state ###
        # ship
        if self.solve_method == "euler":
            ds = MMG.ode_rhs(ship_state, ship_action,ship_wind, update_params)
            ship_state_n = ship_state + self.dt_sim * ds
        elif self.solve_method == "rk4":
            k1 = MMG.ode_rhs(ship_state, ship_action, ship_wind, update_params)
            k2 = MMG.ode_rhs(ship_state + 0.5 * k1 * self.dt_sim, ship_action, ship_wind,update_params)
            k3 = MMG.ode_rhs(ship_state + 0.5 * k2 * self.dt_sim, ship_action, ship_wind,update_params)
            k4 = MMG.ode_rhs(ship_state + 1.0 * k3 * self.dt_sim, ship_action, ship_wind,update_params)
            ds = (1.0 * k1 + 2.0 * k2 + 2.0 * k3 + 1.0 * k4) / 6.0
            ship_state_n = ship_state + self.dt_sim * ds
        else:
            print("solve_method not exist")
            # break
            

        return t_n, ship_state_n

