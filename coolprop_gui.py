# -*- coding: utf-8 -*-
import tkinter as tk
from tkinter import ttk, messagebox
import CoolProp.CoolProp as CP
import sys
import time
import traceback
import webbrowser
from collections import namedtuple
import configparser
import os
import logging
from math import isnan

# --- Pint Import ---
try:
    import pint
except ImportError:
    messagebox.showerror("Dependency Error", "The 'pint' library is required.\nPlease install it using: pip install pint")
    sys.exit(1)

# --- Plotting Imports ---
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import warnings
warnings.filterwarnings("ignore", category=matplotlib.MatplotlibDeprecationWarning)

# --- Constants ---
APP_VERSION = "1.16" # P-T Saturation Added
CONFIG_FILE = "coolprop_gui_config.ini"
COOLPROP_DOCS_URL = "http://www.coolprop.org/coolprop/HighLevelAPI.html#parameter-table"
DEFAULT_UNIT = '?'
ERROR_MSG_PREFIX = "Error: "
WARNING_MSG_PREFIX = "Warning: "
INFO_MSG_PREFIX = "Info: "
FLUID_INFO_PREFIX = "Fluid Information"
PHASE_PREFIX = "Phase: "
COPY_BUTTON_TEXT = "Copy Results"
COPIED_BUTTON_TEXT = "Copied!"
DEFAULT_PHASE_TEXT = PHASE_PREFIX + "-"
NA_TEXT = "N/A"
COOLPROP_ERROR_PREFIX = "CoolProp Error: "
RANKINE_TAB_TEXT = " Rankine Cycle "
REFRIG_TAB_TEXT = " Refrigeration Cycle "
PSYCHRO_TAB_TEXT = " Psychrometry"
PSAT_TAB_TEXT = " P-T Saturation " # NEW Constant

# --- Logging Setup ---
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# --- Pint Unit Registry ---
ureg = pint.UnitRegistry(autoconvert_offset_to_baseunit=True)
Q_ = ureg.Quantity

# --- Pint Unit Options Setup ---
PINT_UNIT_OPTIONS = {
    "Temperature": ['kelvin', 'degC', 'degF'], "Pressure": ['pascal', 'kilopascal', 'bar', 'psi'],
    "Density": ['kg/m**3', 'lb/ft**3'], "SpecificEnergy": ['J/kg', 'kJ/kg', 'BTU/lb'],
    "SpecificEntropy": ['J/(kg*K)', 'kJ/(kg*K)', 'BTU/(lb*degR)'],
    "Conductivity": ['W/(m*K)', 'BTU/(ft*hr*degF)'], "Viscosity": ['Pa*s', 'cP'],
    "SpeedOfSound": ['m/s'], "SurfaceTension": ['N/m'], "MolarMass": ['kg/mol', 'lb/mol'],
    "Dimensionless": ['dimensionless'], "Default": ['dimensionless']
}
PropertyDetails = namedtuple("PropertyDetails", ["prop_type", "si_unit_str", "description"])

PROPERTY_INFO = {
    'T': PropertyDetails('Temperature', 'kelvin', 'Temperature'), 'P': PropertyDetails('Pressure', 'pascal', 'Pressure'),
    'D': PropertyDetails('Density', 'kg/m**3', 'Mass Density'), 'DMASS': PropertyDetails('Density', 'kg/m**3', 'Mass Density'),
    'H': PropertyDetails('SpecificEnergy', 'J/kg', 'Mass Specific Enthalpy'), 'HMASS': PropertyDetails('SpecificEnergy', 'J/kg', 'Mass Specific Enthalpy'),
    'S': PropertyDetails('SpecificEntropy', 'J/(kg*kelvin)', 'Mass Specific Entropy'), 'SMASS': PropertyDetails('SpecificEntropy', 'J/(kg*kelvin)', 'Mass Specific Entropy'),
    'U': PropertyDetails('SpecificEnergy', 'J/kg', 'Mass Specific Internal Energy'), 'UMASS': PropertyDetails('SpecificEnergy', 'J/kg', 'Mass Specific Internal Energy'),
    'Q': PropertyDetails('Dimensionless', 'dimensionless', 'Quality'),
    'CPMASS': PropertyDetails('SpecificEntropy', 'J/(kg*kelvin)', 'Mass Specific Heat (Const P)'),
    'CVMASS': PropertyDetails('SpecificEntropy', 'J/(kg*kelvin)', 'Mass Specific Heat (Const V)'),
    'L': PropertyDetails('Conductivity', 'W/(m*kelvin)', 'Thermal Conductivity'),
    'CONDUCTIVITY': PropertyDetails('Conductivity', 'W/(m*kelvin)', 'Thermal Conductivity'),
    'V': PropertyDetails('Viscosity', 'Pa*s', 'Viscosity'),
    'VISCOSITY': PropertyDetails('Viscosity', 'Pa*s', 'Viscosity'),
    'A': PropertyDetails('SpeedOfSound', 'm/s', 'Speed of Sound'),
    'SPEED_OF_SOUND': PropertyDetails('SpeedOfSound', 'm/s', 'Speed of Sound'),
    'I': PropertyDetails('SurfaceTension', 'N/m', 'Surface Tension'),
    'SURFACE_TENSION': PropertyDetails('SurfaceTension', 'N/m', 'Surface Tension'),
    'M': PropertyDetails('MolarMass', 'kg/mol', 'Molar Mass'),
    'MOLARMASS': PropertyDetails('MolarMass', 'kg/mol', 'Molar Mass'),
    'G': PropertyDetails('SpecificEnergy', 'J/kg', 'Mass Specific Gibbs Energy'),
    'GMASS': PropertyDetails('SpecificEnergy', 'J/kg', 'Mass Specific Gibbs Energy'),
    'GIBBSMASS': PropertyDetails('SpecificEnergy', 'J/kg', 'Mass Specific Gibbs Energy'),
    'DEFAULT': PropertyDetails('Default', 'dimensionless', 'Unknown Property')
}

# --- Thermodynamics Calculation Class ---
class ThermodynamicsCalculator:
    """ Handles CoolProp interactions including dome and line calculations """
    def __init__(self):
        self.ureg = ureg
        self.Q_ = Q_

    def get_fluid_list(self):
        try:
            fluids_str = CP.get_global_param_string("fluids_list")
            fluids = sorted(fluids_str.split(','))
            if not fluids or not fluids_str:
                logging.warning("CoolProp returned empty fluid list string. Using minimal fallback.")
                return ["Water", "Air", "R134a"]
            return fluids
        except Exception as e:
            logging.error(f"Could not fetch full fluid list from CoolProp. Using fallback list. Error: {e}", exc_info=True)
            return sorted(["Water", "Air", "R134a", "CO2", "Propane", "Ammonia", "Nitrogen"])

    def get_coolprop_version(self):
        try: return CP.get_global_param_string('version')
        except Exception as e: logging.error(f"Could not get CoolProp version: {e}"); return NA_TEXT

    def get_fluid_info(self, fluid_name):
        info = {'Tcrit': NA_TEXT, 'Pcrit': NA_TEXT, 'M': NA_TEXT, 'Tmin': NA_TEXT, 'Tmax': NA_TEXT}
        param_map = {'Tcrit': 'Tcrit', 'Pcrit': 'Pcrit', 'M': 'M', 'Tmin': 'Tmin', 'Tmax': 'Tmax'}
        if not fluid_name: return info
        for key, cp_key in param_map.items():
            try: info[key] = float(CP.PropsSI(cp_key, '', 0, '', 0, fluid_name))
            except ValueError as e:
                err_str = str(e).lower()
                if not any(s in err_str for s in ["single", "undefined", "unable", "critical", "triple", "incompressible"]):
                       logging.warning(f"Fluid info ValueError for '{cp_key}' on fluid '{fluid_name}': {e}")
                info[key] = NA_TEXT
            except Exception as e:
                logging.error(f"Fluid info Exception for '{cp_key}' on fluid '{fluid_name}': {e}", exc_info=True)
                info[key] = NA_TEXT
        return info

    def _parse_coolprop_error(self, e, context=""):
        err_str = str(e); prefix = f"{COOLPROP_ERROR_PREFIX}{context}: " if context else f"{COOLPROP_ERROR_PREFIX}"
        if "Input value [" in err_str and "is outside the range" in err_str:
            try: parts = err_str.split("Input value [")[1].split("] is outside the range")[0]; return f"{prefix}Input value {parts} is outside the valid range."
            except IndexError: pass
        if "PropsSI:" in err_str: detail = err_str.split("PropsSI:",1)[-1].strip(); return f"{prefix}{detail}"
        if "PhaseSI:" in err_str: detail = err_str.split("PhaseSI:",1)[-1].strip(); return f"{prefix}{detail}"
        if "Unable to solve" in err_str: return f"{prefix}Unable to solve for the given state."
        if "incompressible" in err_str.lower(): return f"{prefix}Property likely not available for incompressible fluid."
        logging.warning(f"Could not parse CoolProp error: {err_str}")
        return f"{prefix}{err_str}"

    def calculate_properties(self, fluid, p1c, p1v, p2c, p2v, codes_to_calc):
        results = {"error": None, "phase": "Unknown"}; calc_vals = {}; general_error_msg = None
        try:
            phase_str = CP.PhaseSI(p1c, p1v, p2c, p2v, fluid); results["phase"] = phase_str if phase_str else "Unknown"
            logging.info(f"Phase calculated: {results['phase']} for {fluid} with {p1c}={p1v}, {p2c}={p2v}")
        except ValueError as e:
            parsed_err = self._parse_coolprop_error(e, "Phase Calculation"); logging.warning(f"Phase calculation ValueError for {fluid}: {e}")
            results["phase"] = "Phase Error"; general_error_msg = parsed_err; results["error"] = general_error_msg
            results["values"] = {c: general_error_msg for c in codes_to_calc}; return results
        except Exception as e:
            logging.error(f"Unexpected Phase calculation Exception for {fluid}: {e}", exc_info=True); results["phase"] = "Phase Error"
            general_error_msg = f"{ERROR_MSG_PREFIX}Unexpected error during phase calculation."; results["error"] = general_error_msg
            results["values"] = {c: general_error_msg for c in codes_to_calc}; return results
        for code in codes_to_calc:
            if code == p1c: calc_vals[code] = p1v; continue
            if code == p2c: calc_vals[code] = p2v; continue
            try: value = CP.PropsSI(code, p1c, p1v, p2c, p2v, fluid); calc_vals[code] = value
            except ValueError as e:
                parsed_err = self._parse_coolprop_error(e, f"Property '{code}'"); logging.warning(f"PropsSI ValueError for '{code}' on {fluid}: {e}")
                calc_vals[code] = parsed_err;
                if general_error_msg is None: general_error_msg = parsed_err; results["error"] = general_error_msg
            except Exception as e:
                logging.error(f"PropsSI Exception for '{code}' on {fluid}: {e}", exc_info=True); error_msg = f"{ERROR_MSG_PREFIX}Unexpected error calculating '{code}'."
                calc_vals[code] = error_msg
                if general_error_msg is None: general_error_msg = error_msg; results["error"] = general_error_msg
        results["values"] = calc_vals
        if general_error_msg and results["error"] is None: results["error"] = general_error_msg
        if results["error"] and results["phase"] != "Phase Error": logging.warning(f"Calculation completed with errors for {fluid}. First error: {results['error']}")
        return results

    # --- Dome Calculation Methods ---
    def calculate_ph_dome(self, fluid_name):
        dome_data = {'p': None, 'h_liq': None, 'h_vap': None, 'crit_point': None, 'error': None}
        logging.info(f"Calculating P-h dome for {fluid_name}...")
        N = 100
        try:
            pcrit_si = CP.PropsSI('Pcrit', fluid_name)
            tcrit_si = CP.PropsSI('Tcrit', fluid_name)
            try: pmin_si = CP.PropsSI('Pmin', fluid_name)
            except ValueError: pmin_si = CP.PropsSI('ptriple', fluid_name)
            p_start = pmin_si * 1.01; p_end = pcrit_si * 0.999
            if p_start >= p_end or not np.isfinite(p_start) or not np.isfinite(p_end):
                raise ValueError(f"Invalid or non-finite pressure range for dome calculation [{p_start=:.3g} Pa to {p_end=:.3g} Pa]")
            p_range_si = np.logspace(np.log10(p_start), np.log10(p_end), N)
            h_l_si = np.full(N, np.nan); h_v_si = np.full(N, np.nan); valid_indices = []
            for i, p_si_val in enumerate(p_range_si):
                try:
                    h_l_si[i] = CP.PropsSI('H','P', p_si_val,'Q', 0, fluid_name)
                    h_v_si[i] = CP.PropsSI('H','P', p_si_val,'Q', 1, fluid_name)
                    valid_indices.append(i)
                except ValueError: continue
            if not valid_indices: raise ValueError("No valid saturation points calculated for the P-h dome.")
            valid_mask = np.array(valid_indices)
            p_plot_arr = self.Q_(p_range_si[valid_mask], 'pascal').to('kPa').m
            h_l_plot = self.Q_(h_l_si[valid_mask], 'J/kg').to('kJ/kg').m
            h_v_plot = self.Q_(h_v_si[valid_mask], 'J/kg').to('kJ/kg').m
            if len(p_plot_arr) < 2: raise ValueError(f"Insufficient valid points ({len(p_plot_arr)}) for P-h dome line.")
            dome_data['p'] = p_plot_arr; dome_data['h_liq'] = h_l_plot; dome_data['h_vap'] = h_v_plot
            try:
                hcrit_si = CP.PropsSI('H','T', tcrit_si,'P', pcrit_si, fluid_name)
                hc_plot = self.Q_(hcrit_si,'J/kg').to('kJ/kg').m; pc_plot = self.Q_(pcrit_si,'pascal').to('kPa').m
                dome_data['crit_point'] = (hc_plot, pc_plot)
            except ValueError as e_crit:
                logging.warning(f"Could not calculate critical point enthalpy for P-h dome of {fluid_name}: {e_crit}"); dome_data['crit_point'] = None
        except ValueError as e_calc:
             error_msg = f"P-h Dome Calculation Error: {e_calc}"; logging.error(f"{error_msg} for fluid {fluid_name}", exc_info=False); dome_data['error'] = error_msg
        except Exception as e_glob:
             error_msg = f"Unexpected P-h Dome Error: {e_glob}"; logging.error(f"{error_msg} for fluid {fluid_name}", exc_info=True); dome_data['error'] = error_msg
        return dome_data

    def calculate_ts_dome(self, fluid_name):
        dome_data = {'t': None, 's_liq': None, 's_vap': None, 'crit_point': None, 'error': None}
        logging.info(f"Calculating T-s dome for {fluid_name}...")
        N = 100
        try:
            tcrit_si = CP.PropsSI('Tcrit', fluid_name); pcrit_si = CP.PropsSI('Pcrit', fluid_name)
            try: tmin_si = CP.PropsSI('Tmin', fluid_name)
            except ValueError: tmin_si = CP.PropsSI('Ttriple', fluid_name)
            t_start = tmin_si * 1.01; t_end = tcrit_si * 0.999
            if t_start >= t_end or not np.isfinite(t_start) or not np.isfinite(t_end):
                raise ValueError(f"Invalid or non-finite temperature range for dome calculation [{t_start=:.3g} K to {t_end=:.3g} K]")
            t_range_si = np.linspace(t_start, t_end, N)
            s_l_si = np.full(N, np.nan); s_v_si = np.full(N, np.nan); valid_indices = []
            for i, t_si_val in enumerate(t_range_si):
                try:
                    s_l_si[i] = CP.PropsSI('S', 'T', t_si_val, 'Q', 0, fluid_name)
                    s_v_si[i] = CP.PropsSI('S', 'T', t_si_val, 'Q', 1, fluid_name)
                    valid_indices.append(i)
                except ValueError: continue
            if not valid_indices: raise ValueError("No valid saturation points calculated for the T-s dome.")
            valid_mask = np.array(valid_indices)
            t_plot_arr = self.Q_(t_range_si[valid_mask], 'kelvin').to('K').m
            s_l_plot = self.Q_(s_l_si[valid_mask], 'J/(kg*K)').to('kJ/(kg*K)').m
            s_v_plot = self.Q_(s_v_si[valid_mask], 'J/(kg*K)').to('kJ/(kg*K)').m
            if len(t_plot_arr) < 2: raise ValueError(f"Insufficient valid points ({len(t_plot_arr)}) for T-s dome line.")
            dome_data['t'] = t_plot_arr; dome_data['s_liq'] = s_l_plot; dome_data['s_vap'] = s_v_plot
            try:
                scrit_si = CP.PropsSI('S','T', tcrit_si,'P', pcrit_si, fluid_name)
                sc_plot = self.Q_(scrit_si, 'J/(kg*K)').to('kJ/(kg*K)').m; tc_plot = self.Q_(tcrit_si, 'kelvin').to('K').m
                dome_data['crit_point'] = (sc_plot, tc_plot)
            except ValueError as e_crit:
                logging.warning(f"Could not calculate critical point entropy for T-s dome of {fluid_name}: {e_crit}"); dome_data['crit_point'] = None
        except ValueError as e_calc:
             error_msg = f"T-s Dome Calculation Error: {e_calc}"; logging.error(f"{error_msg} for fluid {fluid_name}", exc_info=False); dome_data['error'] = error_msg
        except Exception as e_glob:
             error_msg = f"Unexpected T-s Dome Error: {e_glob}"; logging.error(f"{error_msg} for fluid {fluid_name}", exc_info=True); dome_data['error'] = error_msg
        return dome_data

    # --- Isotherm/Isobar Calculation Methods ---
    def calculate_ph_isotherm(self, fluid_name, temp_k, num_points=50):
        isotherm_data = {'p': [], 'h': [], 'error': None}
        logging.info(f"Calculating P-h isotherm for {fluid_name} at T={temp_k:.2f} K...")
        try:
            pcrit_si = CP.PropsSI('Pcrit', fluid_name); tcrit_si = CP.PropsSI('Tcrit', fluid_name)
            try: pmin_si = CP.PropsSI('Pmin', fluid_name)
            except ValueError: pmin_si = CP.PropsSI('ptriple', fluid_name)
            pmax_plot_factor = 5; pmin_si = max(pmin_si, 1e-3)
            if not np.isfinite(pmin_si): pmin_si = 1.0
            if temp_k >= tcrit_si:
                p_end_si = pcrit_si * pmax_plot_factor
                if pmin_si * 1.01 >= p_end_si: p_range_si = np.linspace(pmin_si, p_end_si, num_points)
                else: p_range_si = np.logspace(np.log10(pmin_si * 1.01), np.log10(p_end_si), num_points)
            else:
                try: psat_si = CP.PropsSI('P', 'T', temp_k, 'Q', 0, fluid_name)
                except ValueError: raise ValueError(f"Cannot determine saturation pressure at T={temp_k:.2f} K")
                p_start_liq = pmin_si * 1.01; p_end_liq = psat_si * 0.999
                if p_start_liq >= p_end_liq: p_range_liq_si = np.array([p_start_liq])
                else: p_range_liq_si = np.logspace(np.log10(p_start_liq), np.log10(p_end_liq), num_points // 2 + (num_points % 2))
                p_start_vap = psat_si * 1.001; p_end_vap = pcrit_si * pmax_plot_factor
                if p_start_vap >= p_end_vap: p_range_vap_si = np.array([p_start_vap])
                else: p_range_vap_si = np.logspace(np.log10(p_start_vap), np.log10(p_end_vap), num_points // 2)
                p_range_si = np.concatenate((p_range_liq_si, p_range_vap_si))
            h_vals_si = []; p_vals_si_valid = []
            for p_si in p_range_si:
                try:
                    if p_si <= 0 or not np.isfinite(p_si): continue
                    h_si = CP.PropsSI('H', 'P', p_si, 'T', temp_k, fluid_name)
                    if np.isfinite(h_si): h_vals_si.append(h_si); p_vals_si_valid.append(p_si)
                except ValueError: continue
            if not p_vals_si_valid: raise ValueError("No valid points calculated for the isotherm.")
            isotherm_data['p'] = self.Q_(np.array(p_vals_si_valid), 'pascal').to('kPa').m
            isotherm_data['h'] = self.Q_(np.array(h_vals_si), 'J/kg').to('kJ/kg').m
        except ValueError as e_calc:
            error_msg = f"P-h Isotherm Calc Error (T={temp_k:.2f} K): {e_calc}"; logging.error(f"{error_msg} for fluid {fluid_name}", exc_info=False); isotherm_data['error'] = error_msg
        except Exception as e_glob:
            error_msg = f"Unexpected P-h Isotherm Error (T={temp_k:.2f} K): {e_glob}"; logging.error(f"{error_msg} for fluid {fluid_name}", exc_info=True); isotherm_data['error'] = error_msg
        return isotherm_data

    def calculate_ts_isobar(self, fluid_name, press_pa, num_points=50):
        isobar_data = {'t': [], 's': [], 'error': None}
        logging.info(f"Calculating T-s isobar for {fluid_name} at P={press_pa/1e3:.2f} kPa...")
        try:
            tcrit_si = CP.PropsSI('Tcrit', fluid_name); pcrit_si = CP.PropsSI('Pcrit', fluid_name)
            try: tmin_si = CP.PropsSI('Tmin', fluid_name)
            except ValueError: tmin_si = CP.PropsSI('Ttriple', fluid_name)
            tmax_si = CP.PropsSI('Tmax', fluid_name)
            if not np.isfinite(tmin_si): tmin_si = 200.0
            if not np.isfinite(tmax_si): tmax_si = 2000.0
            tmin_si = max(tmin_si, 1e-2)
            if press_pa >= pcrit_si:
                t_start = tmin_si * 1.01; t_end = tmax_si * 0.99
                if t_start >= t_end: raise ValueError(f"Invalid temperature range for supercritical isobar [{t_start=:.2f}K, {t_end=:.2f}K]")
                t_range_si = np.linspace(t_start, t_end, num_points)
            else:
                try: tsat_si = CP.PropsSI('T', 'P', press_pa, 'Q', 0, fluid_name)
                except ValueError: raise ValueError(f"Cannot determine saturation temperature at P={press_pa/1e3:.2f} kPa")
                t_start_liq = tmin_si * 1.01; t_end_liq = tsat_si * 0.999
                if t_start_liq >= t_end_liq: t_range_liq_si = np.array([t_start_liq])
                else: t_range_liq_si = np.linspace(t_start_liq, t_end_liq, num_points // 2 + (num_points % 2))
                t_start_vap = tsat_si * 1.001; t_end_vap = tmax_si * 0.99
                if t_start_vap >= t_end_vap: t_range_vap_si = np.array([t_start_vap])
                else: t_range_vap_si = np.linspace(t_start_vap, t_end_vap, num_points // 2)
                t_range_si = np.concatenate((t_range_liq_si, t_range_vap_si))
            s_vals_si = []; t_vals_si_valid = []
            for t_si in t_range_si:
                try:
                    if t_si <= 0 or not np.isfinite(t_si): continue
                    s_si = CP.PropsSI('S', 'T', t_si, 'P', press_pa, fluid_name)
                    if np.isfinite(s_si): s_vals_si.append(s_si); t_vals_si_valid.append(t_si)
                except ValueError: continue
            if not t_vals_si_valid: raise ValueError("No valid points calculated for the isobar.")
            isobar_data['t'] = self.Q_(np.array(t_vals_si_valid), 'kelvin').to('K').m
            isobar_data['s'] = self.Q_(np.array(s_vals_si), 'J/(kg*K)').to('kJ/(kg*K)').m
        except ValueError as e_calc:
            error_msg = f"T-s Isobar Calc Error (P={press_pa/1e3:.2f} kPa): {e_calc}"; logging.error(f"{error_msg} for fluid {fluid_name}", exc_info=False); isobar_data['error'] = error_msg
        except Exception as e_glob:
            error_msg = f"Unexpected T-s Isobar Error (P={press_pa/1e3:.2f} kPa): {e_glob}"; logging.error(f"{error_msg} for fluid {fluid_name}", exc_info=True); isobar_data['error'] = error_msg
        return isobar_data

    # --- Rankine Cycle Calculation ---
    def calculate_rankine_cycle(self, fluid_name, p_boiler_pa, t_boiler_k, p_cond_pa, eta_turbine):
        results = {'states': {}, 'metrics': {}, 'error': None}; state_props = {}
        try:
            if not (0.0 <= eta_turbine <= 1.0): raise ValueError("Turbine efficiency must be between 0.0 and 1.0")
            if p_boiler_pa <= p_cond_pa: raise ValueError("Boiler pressure must be greater than condenser pressure")
            try:
                 t_sat_boiler = CP.PropsSI('T', 'P', p_boiler_pa, 'Q', 1, fluid_name)
                 if t_boiler_k < t_sat_boiler - 0.01: raise ValueError(f"Boiler temperature ({t_boiler_k:.2f} K) is below saturation temperature ({t_sat_boiler:.2f} K) at boiler pressure.")
            except ValueError as e: logging.warning(f"Could not verify boiler saturation temp: {e}")

            logging.info(f"Rankine State 1: P={p_cond_pa:.1f} Pa, Q=0")
            h1 = CP.PropsSI('H', 'P', p_cond_pa, 'Q', 0, fluid_name); s1 = CP.PropsSI('S', 'P', p_cond_pa, 'Q', 0, fluid_name); t1 = CP.PropsSI('T', 'P', p_cond_pa, 'Q', 0, fluid_name)
            state_props[1] = {'T': t1, 'P': p_cond_pa, 'H': h1, 'S': s1, 'Q': 0.0}
            if any(isnan(val) for val in state_props[1].values() if isinstance(val, float)): raise ValueError("Calculation failed at State 1 (Pump Inlet)")

            s2 = s1
            logging.info(f"Rankine State 2: P={p_boiler_pa:.1f} Pa, S={s2:.2f} J/kg/K (ideal pump)")
            h2 = CP.PropsSI('H', 'P', p_boiler_pa, 'S', s2, fluid_name); t2 = CP.PropsSI('T', 'P', p_boiler_pa, 'S', s2, fluid_name)
            state_props[2] = {'T': t2, 'P': p_boiler_pa, 'H': h2, 'S': s2, 'Q': np.nan}
            if any(isnan(val) for val in [t2, h2] if isinstance(val, float)): raise ValueError("Calculation failed at State 2 (Pump Outlet - Ideal)")
            try: state_props[2]['Q'] = CP.PropsSI('Q', 'P', p_boiler_pa, 'H', h2, fluid_name)
            except: pass

            logging.info(f"Rankine State 3: P={p_boiler_pa:.1f} Pa, T={t_boiler_k:.2f} K")
            h3 = CP.PropsSI('H', 'P', p_boiler_pa, 'T', t_boiler_k, fluid_name); s3 = CP.PropsSI('S', 'P', p_boiler_pa, 'T', t_boiler_k, fluid_name)
            state_props[3] = {'T': t_boiler_k, 'P': p_boiler_pa, 'H': h3, 'S': s3, 'Q': np.nan}
            if any(isnan(val) for val in [h3, s3] if isinstance(val, float)): raise ValueError("Calculation failed at State 3 (Turbine Inlet)")
            try: state_props[3]['Q'] = CP.PropsSI('Q', 'P', p_boiler_pa, 'H', h3, fluid_name)
            except: pass

            s4s = s3; logging.info(f"Rankine State 4s: P={p_cond_pa:.1f} Pa, S={s4s:.2f} J/kg/K")
            h4s = CP.PropsSI('H', 'P', p_cond_pa, 'S', s4s, fluid_name)
            if isnan(h4s): raise ValueError("Calculation failed at State 4s (Ideal Turbine Outlet)")
            w_turbine_ideal = h3 - h4s; w_turbine_actual = w_turbine_ideal * eta_turbine; h4 = h3 - w_turbine_actual
            logging.info(f"Rankine State 4: P={p_cond_pa:.1f} Pa, H={h4:.2f} J/kg (Actual, eta={eta_turbine:.2f})")
            t4 = CP.PropsSI('T', 'P', p_cond_pa, 'H', h4, fluid_name); s4 = CP.PropsSI('S', 'P', p_cond_pa, 'H', h4, fluid_name); q4 = CP.PropsSI('Q', 'P', p_cond_pa, 'H', h4, fluid_name)
            state_props[4] = {'T': t4, 'P': p_cond_pa, 'H': h4, 'S': s4, 'Q': q4}
            if any(isnan(val) for val in [t4, s4, q4] if isinstance(val, float)): logging.warning(f"NaN detected in calculated State 4 properties. Proceeding cautiously.")

            w_pump_actual = h2 - h1; q_boiler = h3 - h2; q_condenser = h4 - h1
            if q_boiler <= 0: logging.warning(f"Boiler heat input is zero or negative. Check states."); thermal_efficiency = 0.0; results['error'] = "Warning: Boiler heat input is not positive."
            else: w_net = w_turbine_actual - w_pump_actual; thermal_efficiency = w_net / q_boiler

            results['states'] = state_props
            results['metrics'] = {'W_pump': w_pump_actual, 'Q_boiler': q_boiler, 'W_turbine': w_turbine_actual, 'Q_condenser': q_condenser, 'eta_thermal': thermal_efficiency}
            logging.info(f"Rankine cycle calculated successfully. Efficiency = {thermal_efficiency*100:.2f}%")
        except ValueError as ve:
            error_msg = f"Rankine Calc Error: {ve}"; logging.error(f"{error_msg} for fluid {fluid_name}", exc_info=False); results['error'] = error_msg
        except Exception as e:
            error_msg = f"Unexpected Rankine Calc Error: {e}"; logging.error(f"{error_msg} for fluid {fluid_name}", exc_info=True); results['error'] = error_msg
        return results

    # --- NEW: P-T Saturation Curve Calculation ---
    def calculate_psat_vs_t(self, fluid_name, num_points=100):
        """ Calculates P-T saturation line from triple point to critical point. """
        psat_data = {'t': None, 'p': None, 'error': None}
        logging.info(f"Calculating P-T saturation curve for {fluid_name}...")
        try:
            try: tmin_si = CP.PropsSI('Ttriple', fluid_name)
            except ValueError: logging.warning(f"Could not get Ttriple for {fluid_name}, using Tmin."); tmin_si = CP.PropsSI('Tmin', fluid_name)
            tcrit_si = CP.PropsSI('Tcrit', fluid_name)
            t_start = tmin_si * 1.001; t_end = tcrit_si * 0.999
            if t_start >= t_end or not np.isfinite(t_start) or not np.isfinite(t_end): raise ValueError(f"Invalid temperature range for P-T curve.")
            t_range_si = np.linspace(t_start, t_end, num_points)
            p_sat_si = np.full(num_points, np.nan); valid_indices = []
            for i, t_si_val in enumerate(t_range_si):
                try:
                    p_sat_si[i] = CP.PropsSI('P', 'T', t_si_val, 'Q', 0, fluid_name); valid_indices.append(i)
                except ValueError: continue
            if not valid_indices: raise ValueError("No valid saturation points calculated for the P-T curve.")
            valid_mask = np.array(valid_indices)
            t_plot_arr = self.Q_(t_range_si[valid_mask], 'kelvin').to('K').m
            p_plot_arr = self.Q_(p_sat_si[valid_mask], 'pascal').to('kPa').m
            if len(t_plot_arr) < 2: raise ValueError(f"Insufficient valid points for P-T curve.")
            psat_data['t'] = t_plot_arr; psat_data['p'] = p_plot_arr
            logging.info(f"P-T saturation curve calculated successfully ({len(t_plot_arr)} points).")
        except ValueError as e_calc:
             error_msg = f"P-T Curve Calc Error: {e_calc}"; logging.error(f"{error_msg} for fluid {fluid_name}", exc_info=False); psat_data['error'] = error_msg
        except Exception as e_glob:
             error_msg = f"Unexpected P-T Curve Error: {e_glob}"; logging.error(f"{error_msg} for fluid {fluid_name}", exc_info=True); psat_data['error'] = error_msg
        return psat_data
    # --- NEW: Refrigeration Cycle Calculation ---
    def calculate_refrigeration_cycle(self, fluid_name, p_evap_pa, p_cond_pa, superheat_k, subcooling_k, eta_compressor):
        """
        Calculates a Vapor Compression Refrigeration Cycle.
        States:
        1: Compressor Inlet (Evap P + Superheat)
        2: Compressor Outlet (Cond P, Efficiency applied)
        3: Condenser Outlet (Cond P - Subcooling)
        4: Expansion Outlet (Isenthalpic to Evap P)
        """
        results = {'states': {}, 'metrics': {}, 'error': None}; state_props = {}
        try:
            # --- Validation ---
            if not (0.0 < eta_compressor <= 1.0): raise ValueError("Efficiency must be between 0 and 1")
            if p_evap_pa >= p_cond_pa: raise ValueError("Evaporator pressure must be lower than condenser pressure")
            if superheat_k < 0 or subcooling_k < 0: raise ValueError("Superheat and Subcooling must be non-negative")

            # --- State 1: Compressor Suction (Superheated Vapor) ---
            # T_sat at Evap Pressure
            try: t_sat_evap = CP.PropsSI('T', 'P', p_evap_pa, 'Q', 1, fluid_name)
            except ValueError: raise ValueError(f"Could not calculate saturation T at P_evap={p_evap_pa:.0f} Pa")
            
            t1 = t_sat_evap + superheat_k
            # Check if T1 exceeds max temp? (Optional robustness)
            h1 = CP.PropsSI('H', 'P', p_evap_pa, 'T', t1, fluid_name)
            s1 = CP.PropsSI('S', 'P', p_evap_pa, 'T', t1, fluid_name)
            state_props[1] = {'T': t1, 'P': p_evap_pa, 'H': h1, 'S': s1, 'Q': 1.0} # Superheated
            logging.info(f"Refrig State 1: T={t1:.2f} K, P={p_evap_pa:.0f} Pa")

            # --- State 2: Compressor Discharge (Superheated) ---
            s2s = s1 # Isentropic
            try: h2s = CP.PropsSI('H', 'P', p_cond_pa, 'S', s2s, fluid_name)
            except ValueError: raise ValueError("Failed to calculate isentropic discharge state (2s)")
            
            w_comp_ideal = h2s - h1
            w_comp_actual = w_comp_ideal / eta_compressor
            h2 = h1 + w_comp_actual
            
            t2 = CP.PropsSI('T', 'P', p_cond_pa, 'H', h2, fluid_name)
            s2 = CP.PropsSI('S', 'P', p_cond_pa, 'H', h2, fluid_name)
            state_props[2] = {'T': t2, 'P': p_cond_pa, 'H': h2, 'S': s2, 'Q': 1.0} # Likely Superheated
            logging.info(f"Refrig State 2: T={t2:.2f} K, H={h2:.0f} J/kg")

            # --- State 3: Condenser Outlet (Subcooled Liquid) ---
            try: t_sat_cond = CP.PropsSI('T', 'P', p_cond_pa, 'Q', 0, fluid_name)
            except ValueError: raise ValueError(f"Could not calculate saturation T at P_cond={p_cond_pa:.0f} Pa")

            t3 = t_sat_cond - subcooling_k
            # Ensure T3 is above freezing/triple point roughly
            try: t_min = CP.PropsSI('Tmin', fluid_name)
            except: t_min = 0
            if t3 < t_min: raise ValueError(f"Calculated Condenser Outlet T ({t3:.2f} K) is below fluid minimum.")

            h3 = CP.PropsSI('H', 'P', p_cond_pa, 'T', t3, fluid_name)
            s3 = CP.PropsSI('S', 'P', p_cond_pa, 'T', t3, fluid_name)
            state_props[3] = {'T': t3, 'P': p_cond_pa, 'H': h3, 'S': s3, 'Q': 0.0} # Subcooled
            logging.info(f"Refrig State 3: T={t3:.2f} K")

            # --- State 4: Expansion Device Outlet (Mixture) ---
            h4 = h3 # Isenthalpic
            p4 = p_evap_pa
            # Calculate properties (likely two-phase)
            t4 = CP.PropsSI('T', 'P', p4, 'H', h4, fluid_name)
            s4 = CP.PropsSI('S', 'P', p4, 'H', h4, fluid_name)
            q4 = CP.PropsSI('Q', 'P', p4, 'H', h4, fluid_name)
            state_props[4] = {'T': t4, 'P': p4, 'H': h4, 'S': s4, 'Q': q4}
            logging.info(f"Refrig State 4: Q={q4:.4f}")

            # --- Metrics ---
            q_evap = h1 - h4
            q_cond = h2 - h3
            # w_comp_actual calculated above
            
            cop_cooling = q_evap / w_comp_actual if w_comp_actual > 0 else 0
            cop_heating = q_cond / w_comp_actual if w_comp_actual > 0 else 0

            results['states'] = state_props
            results['metrics'] = {
                'Q_evap': q_evap,
                'W_comp': w_comp_actual,
                'Q_cond': q_cond,
                'COP_cooling': cop_cooling,
                'COP_heating': cop_heating
            }
            logging.info(f"Refrigeration Cycle Calc Success. COPc={cop_cooling:.2f}")

        except ValueError as ve:
            error_msg = f"Refrig Calc Error: {ve}"
            logging.error(f"{error_msg} for {fluid_name}", exc_info=False)
            results['error'] = error_msg
        except Exception as e:
            error_msg = f"Unexpected Refrig Error: {e}"
            logging.error(f"{error_msg} for {fluid_name}", exc_info=True)
            results['error'] = error_msg
        
        return results
    # --- NEW: Psychrometric (Humid Air) Calculation ---
    def calculate_humid_air_properties(self, t_db_k, rh_frac, p_pa):
        """
        Calculates humid air properties using CoolProp's HAPropsSI.
        Inputs:
            t_db_k: Dry Bulb Temperature [K]
            rh_frac: Relative Humidity [0-1]
            p_pa: Atmospheric Pressure [Pa]
        Returns:
            dict: {'properties': {code: val}, 'error': str | None}
        """
        results = {'properties': {}, 'error': None}
        # Map internal codes to CoolProp HA keys and descriptions
        # HAPropsSI inputs: 'T' (Temp), 'P' (Pressure), 'R' (Rel Humidity)
        outputs_map = {
            'Tdp': 'Dew Point Temperature',
            'Twb': 'Wet Bulb Temperature',
            'W': 'Humidity Ratio',
            'H': 'Mixture Enthalpy',
            'V': 'Mixture Volume',
            'M': 'Mixture Viscosity',
            'K': 'Mixture Conductivity'
        }
        
        logging.info(f"Calculating Psychrometry: T={t_db_k:.2f}K, RH={rh_frac:.2f}, P={p_pa:.0f}Pa")

        try:
            # Validate Inputs before calling C++ code
            if not (0 <= rh_frac <= 1): raise ValueError("Relative Humidity must be between 0 and 1 (0-100%)")
            if t_db_k < 150 or t_db_k > 450: raise ValueError("Temperature outside reasonable psychrometric range")

            calc_props = {}
            for key, desc in outputs_map.items():
                try:
                    # Syntax: HAPropsSI(Output, Input1, Val1, Input2, Val2, Input3, Val3)
                    val = CP.HAPropsSI(key, 'T', t_db_k, 'R', rh_frac, 'P', p_pa)
                    calc_props[key] = val
                except Exception as e:
                    logging.warning(f"HAPropsSI failed for {key}: {e}")
                    calc_props[key] = None
            
            results['properties'] = calc_props
            logging.info("Psychrometric calculation successful.")

        except ValueError as ve:
            error_msg = f"Psychro Input Error: {ve}"
            logging.error(error_msg, exc_info=False)
            results['error'] = error_msg
        except Exception as e:
            error_msg = f"Unexpected Psychro Error: {e}"
            logging.error(error_msg, exc_info=True)
            results['error'] = error_msg
            
        return results

# --- ToolTip Class ---
class ToolTip:
    def __init__(self, widget, text):
        self.widget = widget; self.text = text; self.tooltip_window = None
        self.widget.bind("<Enter>", self.show_tooltip); self.widget.bind("<Leave>", self.hide_tooltip); self._id = None; self._after_id = None
    def show_tooltip(self, event=None):
        self.hide_tooltip();
        if self._after_id: self.widget.after_cancel(self._after_id)
        self._after_id = self.widget.after(750, self._display_tooltip)
    def _display_tooltip(self):
        if self.tooltip_window or not self.text: return
        try:
            x, y, _, _ = self.widget.bbox("insert"); x += self.widget.winfo_rootx() + self.widget.winfo_width() // 2; y += self.widget.winfo_rooty() + self.widget.winfo_height() + 5
            self.tooltip_window = tw = tk.Toplevel(self.widget); tw.wm_overrideredirect(True); tw.wm_geometry(f"+{x}+{y}")
            tw.update_idletasks(); screen_width=tw.winfo_screenwidth(); screen_height=tw.winfo_screenheight(); tooltip_width=tw.winfo_width(); x -= tooltip_width // 2
            if x+tooltip_width > screen_width: x=screen_width-tooltip_width-5
            if x<0: x=5
            if y+tw.winfo_height() > screen_height: y=self.widget.winfo_rooty()-tw.winfo_height()-5
            if y<0: y=5
            tw.wm_geometry(f"+{x}+{y}"); label=tk.Label(tw, text=self.text, justify=tk.LEFT, bg="#ffffe0", relief=tk.SOLID, bd=1, font=("tahoma", "8", "normal")); label.pack(ipadx=1)
        except tk.TclError: self.tooltip_window=None
    def hide_tooltip(self, event=None):
        if self._after_id: self.widget.after_cancel(self._after_id); self._after_id=None
        tw = self.tooltip_window; self.tooltip_window = None
        if tw:
            try: tw.destroy()
            except tk.TclError: pass

# --- Main Application Class ---
class CoolPropApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.ureg = ureg; self.Q_ = Q_; self.calculator = ThermodynamicsCalculator()
        self.title(f"CoolProp Interface v{APP_VERSION} (Pint) - State Point & Plotting")
        self.geometry("950x800"); self.minsize(850, 650)
        self.grid_rowconfigure(0, weight=1); self.grid_columnconfigure(0, weight=1)
        self._apply_style(); self._full_fluid_list = self.calculator.get_fluid_list()
        all_desc = [d.description for c,d in PROPERTY_INFO.items() if c != 'DEFAULT']
        self.property_display_list = sorted(list(set(all_desc)))
        self.display_to_code_map = {}
        for code, details in PROPERTY_INFO.items():
            if code != 'DEFAULT':
                description = details.description
                if description not in self.display_to_code_map or len(code) < len(self.display_to_code_map[description]):
                    self.display_to_code_map[description] = code
        self.code_to_display_map = {v: k for k, v in self.display_to_code_map.items()}
        self.last_si_results = None; self.current_fluid_for_plot = None
        self._setup_variables(); self._create_menu(); self._setup_ui_layout(); self._bind_events(); self._initial_state_setup(); self._load_settings(); self.protocol("WM_DELETE_WINDOW", self._on_closing)

    def _apply_style(self):
        style = ttk.Style(self); available_themes = style.theme_names(); preferred_themes = ["clam", "alt", "default", "vista", "xpnative"]; theme_set=False
        for theme in preferred_themes:
            if theme in available_themes:
                try: style.theme_use(theme); theme_set=True; break;
                except tk.TclError: continue
        if not theme_set:
            try: style.theme_use(available_themes[0]);
            except: pass
        style.map("Treeview", background=[('selected', '#0078D7')], foreground=[('selected', 'white')]); style.configure("Treeview", rowheight=20); style.configure("Treeview.Heading", font=('TkDefaultFont', 9, 'bold')); style.configure("Error.Treeview", foreground='red'); style.map("Error.Treeview", foreground=[('!selected', 'red'), ('selected', 'white')], background=[('selected', '#0078D7')])

    def _setup_variables(self):
        self.fluid_var = tk.StringVar(value="Water"); self.fluid_filter_var = tk.StringVar()
        default_prop1_display = self.code_to_display_map.get("T", "Temperature"); default_prop2_display = self.code_to_display_map.get("P", "Pressure")
        self.prop1_name_var = tk.StringVar(value=default_prop1_display); self.prop1_value_var = tk.StringVar(value="25"); self.prop1_unit_var = tk.StringVar()
        self.prop2_name_var = tk.StringVar(value=default_prop2_display); self.prop2_value_var = tk.StringVar(value="1.01325"); self.prop2_unit_var = tk.StringVar()
        self.phase_var = tk.StringVar(value=DEFAULT_PHASE_TEXT); self.fluid_info_tcrit_var = tk.StringVar(value=NA_TEXT); self.fluid_info_pcrit_var = tk.StringVar(value=NA_TEXT); self.fluid_info_molmass_var = tk.StringVar(value=NA_TEXT); self.fluid_info_tmin_var = tk.StringVar(value=NA_TEXT); self.fluid_info_tmax_var = tk.StringVar(value=NA_TEXT)
        self.show_isotherms_var = tk.BooleanVar(value=False); self.show_isobars_var = tk.BooleanVar(value=False)
        self.rankine_p_boiler_var = tk.StringVar(value="8000"); self.rankine_t_boiler_var = tk.StringVar(value="500"); self.rankine_p_cond_var = tk.StringVar(value="10"); self.rankine_eta_turbine_var = tk.StringVar(value="85")
        self.refrig_p_evap_var = tk.StringVar(value="200") # kPa
        self.refrig_p_cond_var = tk.StringVar(value="1200") # kPa
        self.refrig_sh_var = tk.StringVar(value="5") # K
        self.refrig_sc_var = tk.StringVar(value="5") # K
        self.refrig_eta_var = tk.StringVar(value="80")
        self.psychro_tdb_var = tk.StringVar(value="25") # degC
        self.psychro_rh_var = tk.StringVar(value="50") # %
        self.psychro_p_var = tk.StringVar(value="101.325") # kPa

    def _create_menu(self):
        self.menu_bar=tk.Menu(self); self.config(menu=self.menu_bar); f=tk.Menu(self.menu_bar, tearoff=0); self.menu_bar.add_cascade(label="File", menu=f); f.add_command(label="Exit", command=self._on_closing)
        h=tk.Menu(self.menu_bar, tearoff=0); self.menu_bar.add_cascade(label="Help", menu=h); h.add_command(label="CoolProp Params...", command=self._open_coolprop_docs); h.add_separator(); h.add_command(label="About...", command=self._show_about_dialog)

    def _open_coolprop_docs(self):
        try: webbrowser.open_new_tab(COOLPROP_DOCS_URL)
        except Exception as e: logging.error(f"Could not open link {COOLPROP_DOCS_URL}: {e}", exc_info=True); messagebox.showerror("Error", f"Could not open link:\n{e}", parent=self)

    def _show_about_dialog(self):
        cp_v = self.calculator.get_coolprop_version(); msg=(f"CoolProp GUI v{APP_VERSION}\n\nCoolProp: {cp_v}\nPython: {sys.version.split()[0]}\nTk: {tk.TkVersion}"); messagebox.showinfo("About", msg, parent=self)

    def _setup_ui_layout(self):
        self.notebook = ttk.Notebook(self); self.notebook.grid(row=0, column=0, sticky="nsew", padx=5, pady=5); self.notebook.enable_traversal()
        self.calculator_tab = ttk.Frame(self.notebook, padding="5"); self.ph_plot_tab = ttk.Frame(self.notebook, padding="5"); self.ts_plot_tab = ttk.Frame(self.notebook, padding="5")
        self.notebook.add(self.calculator_tab, text=" Calculator "); self.notebook.add(self.ph_plot_tab, text=" P-h Diagram "); self.notebook.add(self.ts_plot_tab, text=" T-s Diagram ")

        self.rankine_tab = ttk.Frame(self.notebook, padding="5"); self.notebook.add(self.rankine_tab, text=RANKINE_TAB_TEXT)
        self.refrig_tab = ttk.Frame(self.notebook, padding="5")
        self.notebook.add(self.refrig_tab, text=REFRIG_TAB_TEXT)
        self.psychro_tab = ttk.Frame(self.notebook, padding="5")
        self.notebook.add(self.psychro_tab, text=PSYCHRO_TAB_TEXT)
        # --- NEW P-T Saturation Tab Setup ---
        self.psat_plot_tab = ttk.Frame(self.notebook, padding="5")
        self.notebook.add(self.psat_plot_tab, text=PSAT_TAB_TEXT)

        self.rankine_tab.grid_rowconfigure(1, weight=1); self.rankine_tab.grid_columnconfigure(0, weight=1); self.rankine_tab.grid_columnconfigure(1, weight=2)
        self.refrig_tab.grid_rowconfigure(1, weight=1); self.refrig_tab.grid_columnconfigure(0, weight=1); self.refrig_tab.grid_columnconfigure(1, weight=2)
        self.psychro_tab.grid_rowconfigure(0, weight=1); self.psychro_tab.grid_columnconfigure(0, weight=1); self.psychro_tab.grid_columnconfigure(1, weight=2)
        self.calculator_tab.grid_rowconfigure(0, weight=1); self.calculator_tab.grid_columnconfigure(1, weight=2); self.calculator_tab.grid_columnconfigure(0, weight=1, minsize=350)
        self.ph_plot_tab.grid_rowconfigure(0, weight=1); self.ph_plot_tab.grid_columnconfigure(0, weight=1)
        self.ts_plot_tab.grid_rowconfigure(0, weight=1); self.ts_plot_tab.grid_columnconfigure(0, weight=1)
        self.psat_plot_tab.grid_rowconfigure(0, weight=1); self.psat_plot_tab.grid_columnconfigure(0, weight=1)

        self._create_calculator_widgets(self.calculator_tab)
        self._create_plot_widgets()
        self._create_rankine_widgets(self.rankine_tab)
        self._create_refrig_widgets(self.refrig_tab)
        self._create_psychro_widgets(self.psychro_tab)
        self._create_psat_plot_widgets(self.psat_plot_tab) # NEW Call
        self._setup_plot_layout()

    def _create_calculator_widgets(self, parent_frame):
        left_frame = ttk.Frame(parent_frame); left_frame.grid(row=0, column=0, sticky="nsew", padx=(0, 5)); left_frame.grid_rowconfigure(2, weight=1); left_frame.grid_columnconfigure(0, weight=1)
        fluid_frame = ttk.LabelFrame(left_frame, text="Fluid Selection", padding="5"); fluid_frame.grid(row=0, column=0, sticky="new", pady=(0, 5)); fluid_frame.grid_columnconfigure(1, weight=1)
        ttk.Label(fluid_frame, text="Filter:").grid(row=0, column=0, padx=(5,2), pady=3, sticky="w"); self.fluid_filter_entry = ttk.Entry(fluid_frame, textvariable=self.fluid_filter_var, width=20); self.fluid_filter_entry.grid(row=0, column=1, padx=(0,5), pady=3, sticky="ew")
        fluid_list_frame = ttk.Frame(fluid_frame); fluid_list_frame.grid(row=1, column=0, columnspan=2, sticky="nsew", pady=(3,0)); fluid_list_frame.grid_rowconfigure(0, weight=1); fluid_list_frame.grid_columnconfigure(0, weight=1)
        self.fluid_listbox = tk.Listbox(fluid_list_frame, height=6, exportselection=False); fluid_scrollbar = ttk.Scrollbar(fluid_list_frame, orient=tk.VERTICAL, command=self.fluid_listbox.yview); self.fluid_listbox.configure(yscrollcommand=fluid_scrollbar.set); self.fluid_listbox.grid(row=0, column=0, sticky="nsew"); fluid_scrollbar.grid(row=0, column=1, sticky="ns")
        self.input_frame = ttk.LabelFrame(left_frame, text="Input Properties", padding="10"); self.input_frame.grid(row=1, column=0, sticky="new", pady=(5, 5))
        ttk.Label(self.input_frame, text="Prop 1:").grid(row=0, column=0, padx=5, pady=3, sticky="w"); self.prop1_name_combo=ttk.Combobox(self.input_frame, textvariable=self.prop1_name_var, values=self.property_display_list, width=18, state="readonly"); self.prop1_name_combo.grid(row=0, column=1, padx=5, pady=3, sticky="ew"); self.prop1_value_entry=ttk.Entry(self.input_frame, textvariable=self.prop1_value_var, width=12); self.prop1_value_entry.grid(row=0, column=2, padx=5, pady=3, sticky="ew"); self.prop1_unit_combo=ttk.Combobox(self.input_frame, textvariable=self.prop1_unit_var, width=10, state="readonly"); self.prop1_unit_combo.grid(row=0, column=3, padx=(0,5), pady=3, sticky="w")
        ttk.Label(self.input_frame, text="Prop 2:").grid(row=1, column=0, padx=5, pady=3, sticky="w"); self.prop2_name_combo=ttk.Combobox(self.input_frame, textvariable=self.prop2_name_var, values=self.property_display_list, width=18, state="readonly"); self.prop2_name_combo.grid(row=1, column=1, padx=5, pady=3, sticky="ew"); self.prop2_value_entry=ttk.Entry(self.input_frame, textvariable=self.prop2_value_var, width=12); self.prop2_value_entry.grid(row=1, column=2, padx=5, pady=3, sticky="ew"); self.prop2_unit_combo=ttk.Combobox(self.input_frame, textvariable=self.prop2_unit_var, width=10, state="readonly"); self.prop2_unit_combo.grid(row=1, column=3, padx=(0,5), pady=3, sticky="w")
        self.input_frame.grid_columnconfigure(1, weight=1); self.input_frame.grid_columnconfigure(2, weight=1)
        self.output_select_frame=ttk.LabelFrame(left_frame, text="Select Output Properties", padding="10"); self.output_select_frame.grid(row=2, column=0, sticky="nsew", pady=(5,0))
        self.output_listbox=tk.Listbox(self.output_select_frame, selectmode=tk.EXTENDED, height=10, exportselection=False); output_scrollbar=ttk.Scrollbar(self.output_select_frame, orient=tk.VERTICAL, command=self.output_listbox.yview); self.output_listbox.config(yscrollcommand=output_scrollbar.set)
        self.output_listbox.grid(row=0, column=0, sticky="nsew", pady=5); output_scrollbar.grid(row=0, column=1, sticky="ns", pady=5); self.output_select_frame.grid_rowconfigure(0, weight=1); self.output_select_frame.grid_columnconfigure(0, weight=1)
        right_frame = ttk.Frame(parent_frame); right_frame.grid(row=0, column=1, sticky="nsew", padx=(5, 0)); right_frame.grid_rowconfigure(2, weight=1); right_frame.grid_columnconfigure(0, weight=1)
        self.fluid_info_frame=ttk.LabelFrame(right_frame, text="Fluid Information (SI Units)", padding="10"); self.fluid_info_frame.grid(row=0, column=0, sticky="new", pady=(0,5))
        ttk.Label(self.fluid_info_frame, text="Tcrit:").grid(row=0, column=0, padx=5, pady=2, sticky="w"); self.fluid_info_tcrit_label=ttk.Label(self.fluid_info_frame, textvariable=self.fluid_info_tcrit_var, anchor="w"); self.fluid_info_tcrit_label.grid(row=0, column=1, padx=5, pady=2, sticky="ew")
        ttk.Label(self.fluid_info_frame, text="Pcrit:").grid(row=1, column=0, padx=5, pady=2, sticky="w"); self.fluid_info_pcrit_label=ttk.Label(self.fluid_info_frame, textvariable=self.fluid_info_pcrit_var, anchor="w"); self.fluid_info_pcrit_label.grid(row=1, column=1, padx=5, pady=2, sticky="ew")
        ttk.Label(self.fluid_info_frame, text="M:").grid(row=2, column=0, padx=5, pady=2, sticky="w"); self.fluid_info_molmass_label=ttk.Label(self.fluid_info_frame, textvariable=self.fluid_info_molmass_var, anchor="w"); self.fluid_info_molmass_label.grid(row=2, column=1, padx=5, pady=2, sticky="ew")
        ttk.Label(self.fluid_info_frame, text="Tmin:").grid(row=0, column=2, padx=15, pady=2, sticky="w"); self.fluid_info_tmin_label=ttk.Label(self.fluid_info_frame, textvariable=self.fluid_info_tmin_var, anchor="w"); self.fluid_info_tmin_label.grid(row=0, column=3, padx=5, pady=2, sticky="ew")
        ttk.Label(self.fluid_info_frame, text="Tmax:").grid(row=1, column=2, padx=15, pady=2, sticky="w"); self.fluid_info_tmax_label=ttk.Label(self.fluid_info_frame, textvariable=self.fluid_info_tmax_var, anchor="w"); self.fluid_info_tmax_label.grid(row=1, column=3, padx=5, pady=2, sticky="ew")
        self.fluid_info_frame.grid_columnconfigure(1, weight=1); self.fluid_info_frame.grid_columnconfigure(3, weight=1)
        self.button_phase_frame=ttk.Frame(right_frame); self.button_phase_frame.grid(row=1, column=0, sticky="ew", pady=(5,5))
        self.calculate_button=ttk.Button(self.button_phase_frame, text="Calculate", command=self._calculate_property_orchestrator, width=12); self.clear_button=ttk.Button(self.button_phase_frame, text="Clear", command=self._clear_values, width=10); self.phase_label=ttk.Label(self.button_phase_frame, textvariable=self.phase_var, foreground="grey", font="-size 9 -slant italic")
        self.calculate_button.pack(side=tk.LEFT, padx=(0,10), pady=5); self.clear_button.pack(side=tk.LEFT, padx=10, pady=5); self.phase_label.pack(side=tk.LEFT, padx=10, pady=5)
        self.result_frame=ttk.LabelFrame(right_frame, text="Results", padding="10"); self.result_frame.grid(row=2, column=0, sticky="nsew", pady=(5,0)); self.result_frame.grid_rowconfigure(0, weight=1); self.result_frame.grid_columnconfigure(0, weight=1)
        cols = ('prop', 'val', 'unit'); self.result_tree=ttk.Treeview(self.result_frame, columns=cols, show='headings', height=10); self.result_tree.heading('prop', text='Property'); self.result_tree.heading('val', text='Value'); self.result_tree.heading('unit', text='Unit (SI)')
        self.result_tree.column('prop', width=180, anchor=tk.W, stretch=tk.YES); self.result_tree.column('val', width=120, anchor=tk.E, stretch=tk.YES); self.result_tree.column('unit', width=100, anchor=tk.W, stretch=tk.NO)
        result_scroll=ttk.Scrollbar(self.result_frame, orient=tk.VERTICAL, command=self.result_tree.yview); self.result_tree.configure(yscrollcommand=result_scroll.set); self.result_tree.grid(row=0, column=0, sticky='nsew', pady=(0,5)); result_scroll.grid(row=0, column=1, sticky='ns', pady=(0,5))
        self.copy_button=ttk.Button(self.result_frame, text=COPY_BUTTON_TEXT, width=12, command=self._copy_result_to_clipboard); self.copy_button.grid(row=1, column=0, columnspan=2, pady=5, sticky='e')

    def _create_plot_widgets(self):
        self.ph_fig = Figure(figsize=(7, 5), dpi=100); self.ph_ax = self.ph_fig.add_subplot(111)
        self.ph_canvas = FigureCanvasTkAgg(self.ph_fig, master=self.ph_plot_tab)
        self.ph_canvas_widget = self.ph_canvas.get_tk_widget()
        self.ph_toolbar_frame = ttk.Frame(self.ph_plot_tab)
        self.ph_toolbar = NavigationToolbar2Tk(self.ph_canvas, self.ph_toolbar_frame); self.ph_toolbar.update()
        self.isotherms_check = ttk.Checkbutton(self.ph_toolbar_frame, text="Afficher Isotherme(s)", variable=self.show_isotherms_var, command=self._trigger_ph_plot_update)
        self.isotherms_check.pack(side=tk.RIGHT, padx=5); self.ph_toolbar.pack(side=tk.LEFT, fill=tk.X, expand=True)
        self.ph_ax.set_xlabel("Specific Enthalpy (h) [kJ/kg]"); self.ph_ax.set_ylabel("Pressure (P) [kPa]"); self.ph_ax.set_yscale('log'); self.ph_ax.set_title("P-h Diagram"); self.ph_fig.tight_layout()

        self.ts_fig = Figure(figsize=(7, 5), dpi=100); self.ts_ax = self.ts_fig.add_subplot(111)
        self.ts_canvas = FigureCanvasTkAgg(self.ts_fig, master=self.ts_plot_tab)
        self.ts_canvas_widget = self.ts_canvas.get_tk_widget()
        self.ts_toolbar_frame = ttk.Frame(self.ts_plot_tab)
        self.ts_toolbar = NavigationToolbar2Tk(self.ts_canvas, self.ts_toolbar_frame); self.ts_toolbar.update()
        self.isobars_check = ttk.Checkbutton(self.ts_toolbar_frame, text="Afficher Isobare(s)", variable=self.show_isobars_var, command=self._trigger_ts_plot_update)
        self.isobars_check.pack(side=tk.RIGHT, padx=5); self.ts_toolbar.pack(side=tk.LEFT, fill=tk.X, expand=True)
        self.ts_ax.set_xlabel("Entropy [kJ/(kg*K)]"); self.ts_ax.set_ylabel("Temperature [K]"); self.ts_ax.set_title("T-s Diagram"); self.ts_fig.tight_layout()

    # --- NEW: P-T Saturation Plot Widget Creation ---
    def _create_psat_plot_widgets(self, parent_frame):
        self.psat_fig = Figure(figsize=(7, 5), dpi=100)
        self.psat_ax = self.psat_fig.add_subplot(111)
        self.psat_canvas = FigureCanvasTkAgg(self.psat_fig, master=parent_frame)
        self.psat_canvas_widget = self.psat_canvas.get_tk_widget()
        self.psat_toolbar_frame = ttk.Frame(parent_frame)
        self.psat_toolbar = NavigationToolbar2Tk(self.psat_canvas, self.psat_toolbar_frame)
        self.psat_toolbar.update()
        self.psat_toolbar.pack(side=tk.LEFT, fill=tk.X, expand=True)
        self.psat_ax.set_xlabel("Temperature (T) [K]"); self.psat_ax.set_ylabel("Saturation Pressure (P) [kPa]")
        self.psat_ax.set_yscale('log'); self.psat_ax.set_title("P-T Saturation Curve"); self.psat_fig.tight_layout()

    def _create_rankine_widgets(self, parent_frame):
        input_frame = ttk.LabelFrame(parent_frame, text="Rankine Cycle Inputs", padding="10"); input_frame.grid(row=0, column=0, padx=5, pady=5, sticky="new")
        ttk.Label(input_frame, text="Boiler Pressure (P3):").grid(row=0, column=0, padx=5, pady=3, sticky="w"); ttk.Entry(input_frame, textvariable=self.rankine_p_boiler_var, width=10).grid(row=0, column=1, padx=5, pady=3, sticky="ew"); ttk.Label(input_frame, text="kPa").grid(row=0, column=2, padx=(0,5), pady=3, sticky="w")
        ttk.Label(input_frame, text="Boiler Temperature (T3):").grid(row=1, column=0, padx=5, pady=3, sticky="w"); ttk.Entry(input_frame, textvariable=self.rankine_t_boiler_var, width=10).grid(row=1, column=1, padx=5, pady=3, sticky="ew"); ttk.Label(input_frame, text="°C").grid(row=1, column=2, padx=(0,5), pady=3, sticky="w")
        ttk.Label(input_frame, text="Condenser Pressure (P1/P4):").grid(row=2, column=0, padx=5, pady=3, sticky="w"); ttk.Entry(input_frame, textvariable=self.rankine_p_cond_var, width=10).grid(row=2, column=1, padx=5, pady=3, sticky="ew"); ttk.Label(input_frame, text="kPa").grid(row=2, column=2, padx=(0,5), pady=3, sticky="w")
        ttk.Label(input_frame, text="Turbine Efficiency (ηt):").grid(row=3, column=0, padx=5, pady=3, sticky="w"); ttk.Entry(input_frame, textvariable=self.rankine_eta_turbine_var, width=10).grid(row=3, column=1, padx=5, pady=3, sticky="ew"); ttk.Label(input_frame, text="%").grid(row=3, column=2, padx=(0,5), pady=3, sticky="w")
        ToolTip(input_frame.grid_slaves(row=3, column=1)[0], "Isentropic efficiency (0-100)")
        calc_button = ttk.Button(input_frame, text="Calculate Rankine Cycle", command=self._calculate_rankine_orchestrator); calc_button.grid(row=5, column=0, columnspan=3, pady=10)
        results_frame = ttk.LabelFrame(parent_frame, text="Rankine Cycle Results", padding="10"); results_frame.grid(row=0, column=1, rowspan=2, padx=5, pady=5, sticky="nsew"); results_frame.grid_rowconfigure(1, weight=1); results_frame.grid_columnconfigure(0, weight=1)
        metrics_frame = ttk.Frame(results_frame, padding=(0, 5, 0, 10)); metrics_frame.grid(row=0, column=0, sticky="new")
        self.rankine_wpump_label = ttk.Label(metrics_frame, text="W_pump: -"); self.rankine_wpump_label.grid(row=0, column=0, padx=5, pady=1, sticky="w")
        self.rankine_qboiler_label = ttk.Label(metrics_frame, text="Q_boiler: -"); self.rankine_qboiler_label.grid(row=0, column=1, padx=5, pady=1, sticky="w")
        self.rankine_wturbine_label = ttk.Label(metrics_frame, text="W_turbine: -"); self.rankine_wturbine_label.grid(row=1, column=0, padx=5, pady=1, sticky="w")
        self.rankine_qcond_label = ttk.Label(metrics_frame, text="Q_condenser: -"); self.rankine_qcond_label.grid(row=1, column=1, padx=5, pady=1, sticky="w")
        self.rankine_eta_label = ttk.Label(metrics_frame, text="η_thermal: -", font="-weight bold"); self.rankine_eta_label.grid(row=2, column=0, columnspan=2, padx=5, pady=5, sticky="w")
        cols = ('point', 'T', 'P', 'H', 'S', 'Q'); self.rankine_tree = ttk.Treeview(results_frame, columns=cols, show='headings', height=5); self.rankine_tree.grid(row=1, column=0, sticky="nsew")
        self.rankine_tree.heading('point', text='State'); self.rankine_tree.heading('T', text='T [K]'); self.rankine_tree.heading('P', text='P [kPa]'); self.rankine_tree.heading('H', text='H [kJ/kg]'); self.rankine_tree.heading('S', text='S [kJ/kg·K]'); self.rankine_tree.heading('Q', text='Quality')
        self.rankine_tree.column('point', width=50, anchor=tk.CENTER, stretch=tk.NO); self.rankine_tree.column('T', width=80, anchor=tk.E, stretch=tk.YES); self.rankine_tree.column('P', width=90, anchor=tk.E, stretch=tk.YES); self.rankine_tree.column('H', width=100, anchor=tk.E, stretch=tk.YES); self.rankine_tree.column('S', width=100, anchor=tk.E, stretch=tk.YES); self.rankine_tree.column('Q', width=70, anchor=tk.E, stretch=tk.YES)
        scroll = ttk.Scrollbar(results_frame, orient=tk.VERTICAL, command=self.rankine_tree.yview); self.rankine_tree.configure(yscrollcommand=scroll.set); scroll.grid(row=1, column=1, sticky='ns')

    def _setup_plot_layout(self):
        self.ph_canvas_widget.grid(row=0, column=0, sticky="nsew", padx=5, pady=5)
        self.ph_toolbar_frame.grid(row=1, column=0, sticky="ew", padx=5, pady=(0,5))
        self.ts_canvas_widget.grid(row=0, column=0, sticky="nsew", padx=5, pady=5)
        self.ts_toolbar_frame.grid(row=1, column=0, sticky="ew", padx=5, pady=(0,5))
        self.psat_canvas_widget.grid(row=0, column=0, sticky="nsew", padx=5, pady=5)
        self.psat_toolbar_frame.grid(row=1, column=0, sticky="ew", padx=5, pady=(0,5))

    def _bind_events(self):
        self.fluid_filter_entry.bind("<KeyRelease>", self._filter_fluid_list); self.fluid_listbox.bind("<Double-Button-1>", self._on_fluid_select_from_list)
        self.prop1_name_combo.bind("<<ComboboxSelected>>", self._update_unit_options); self.prop2_name_combo.bind("<<ComboboxSelected>>", self._update_unit_options)
        self.prop1_value_entry.bind("<Return>", lambda e: self._calculate_property_orchestrator()); self.prop2_value_entry.bind("<Return>", lambda e: self._calculate_property_orchestrator())
        self.calculate_button.bind("<Return>", lambda e: self._calculate_property_orchestrator()); self.output_listbox.bind("<Double-Button-1>", lambda e: self._calculate_property_orchestrator())
    def _initial_state_setup(self):
        logging.info(f"Initializing CoolPropApp v{APP_VERSION}")
        self._update_unit_options(event=None)
        self._filter_fluid_list()
        default_fluid = self.fluid_var.get()
        if default_fluid:
            try:
                list_items = self.fluid_listbox.get(0, tk.END)
                if default_fluid in list_items:
                    idx = list_items.index(default_fluid)
                    self.fluid_listbox.selection_clear(0, tk.END); self.fluid_listbox.selection_set(idx); self.fluid_listbox.see(idx); self.fluid_filter_var.set(default_fluid)
                    logging.info(f"Default fluid '{default_fluid}' selected in listbox.")
                else: logging.warning(f"Default fluid '{default_fluid}' not found in listbox items."); self.fluid_var.set("")
            except Exception as e: logging.error(f"Error selecting default fluid '{default_fluid}' in listbox: {e}"); self.fluid_var.set("")
        self.output_listbox.delete(0, tk.END)
        for name in self.property_display_list: self.output_listbox.insert(tk.END, name)
        self._reset_output_display()
        self._clear_ph_plot(); self._clear_ts_plot()
        self._clear_psat_plot() # NEW: Clear P-T plot
        self._clear_rankine_results(); self.last_si_results = None; self.current_fluid_for_plot = None
        self._update_fluid_info()
        self._update_psat_plot(self.fluid_var.get()) # NEW: Draw plot for default fluid
        self._create_tooltips()

    def _create_tooltips(self):
        ToolTip(self.fluid_filter_entry, "Type here to filter fluid list below"); ToolTip(self.fluid_listbox, "Double-click to select a fluid from the filtered list")
        ToolTip(self.prop1_name_combo, "Select the first known property"); ToolTip(self.prop1_value_entry, "Enter value for Property 1"); ToolTip(self.prop1_unit_combo, "Select unit for Property 1")
        ToolTip(self.prop2_name_combo, "Select the second known property"); ToolTip(self.prop2_value_entry, "Enter value for Property 2"); ToolTip(self.prop2_unit_combo, "Select unit for Property 2")
        ToolTip(self.output_listbox, "Select outputs (Ctrl+Click, Shift+Click, Dbl-Click to Calculate)"); ToolTip(self.calculate_button, "Calculate selected properties (Enter)"); ToolTip(self.clear_button, "Clear inputs and results"); ToolTip(self.copy_button, "Copy results table")
        ToolTip(self.result_tree, "Calculated properties (SI Units)"); ToolTip(self.phase_label, "Calculated phase"); ToolTip(self.fluid_info_tcrit_label, "Critical temperature (SI units)")
        ToolTip(self.isotherms_check, "Toggle display of isotherm line(s) on P-h diagram"); ToolTip(self.isobars_check, "Toggle display of isobar line(s) on T-s diagram")

    def _load_settings(self):
        logging.info(f"Attempting to load settings from {CONFIG_FILE}..."); config = configparser.ConfigParser()
        if not os.path.exists(CONFIG_FILE): logging.warning(f"Config file '{CONFIG_FILE}' not found. Using defaults."); return
        try:
            config.read(CONFIG_FILE); geometry = config.get('Window', 'Geometry', fallback=None)
            if geometry:
                try: self.geometry(geometry); logging.info(f"  Loaded geometry: {geometry}")
                except tk.TclError as e: logging.warning(f"  Could not apply saved geometry '{geometry}': {e}")
            last_fluid = config.get('State', 'LastFluid', fallback=None)
            if last_fluid and last_fluid in self._full_fluid_list:
                self.fluid_var.set(last_fluid); logging.info(f"  Loaded fluid: {last_fluid}")
                try:
                    list_items = self.fluid_listbox.get(0, tk.END)
                    if last_fluid in list_items:
                        idx = list_items.index(last_fluid); self.fluid_listbox.selection_clear(0, tk.END); self.fluid_listbox.selection_set(idx); self.fluid_listbox.see(idx); self.fluid_filter_var.set(last_fluid); self._update_fluid_info()
                    else: logging.warning(f"  Loaded fluid '{last_fluid}' not found in current listbox items.")
                except Exception as e: logging.warning(f"  Could not select loaded fluid '{last_fluid}' in listbox: {e}")
            elif last_fluid: logging.warning(f"  Saved fluid '{last_fluid}' not found. Keeping default.")
            p1_name = config.get('Inputs', 'Prop1Name', fallback=None); p1_val = config.get('Inputs', 'Prop1Value', fallback=None); p1_unit = config.get('Inputs', 'Prop1Unit', fallback=None)
            p2_name = config.get('Inputs', 'Prop2Name', fallback=None); p2_val = config.get('Inputs', 'Prop2Value', fallback=None); p2_unit = config.get('Inputs', 'Prop2Unit', fallback=None)
            if p1_name and p1_name in self.display_to_code_map:
                self.prop1_name_var.set(p1_name); self._update_unit_options(event=None)
                if p1_unit and p1_unit in self.prop1_unit_combo['values']: self.prop1_unit_var.set(p1_unit)
                if p1_val is not None: self.prop1_value_var.set(p1_val)
            if p2_name and p2_name in self.display_to_code_map:
                self.prop2_name_var.set(p2_name); self._update_unit_options(event=None)
                if p2_unit and p2_unit in self.prop2_unit_combo['values']: self.prop2_unit_var.set(p2_unit)
                if p2_val is not None: self.prop2_value_var.set(p2_val)
            try:
                p1_unit_final = self.prop1_unit_var.get(); p2_unit_final = self.prop2_unit_var.get()
                if p1_unit_final in self.prop1_unit_combo['values']: self.prop1_unit_combo.current(list(self.prop1_unit_combo['values']).index(p1_unit_final))
                if p2_unit_final in self.prop2_unit_combo['values']: self.prop2_unit_combo.current(list(self.prop2_unit_combo['values']).index(p2_unit_final))
            except tk.TclError: pass
            selected_outputs_str = config.get('Outputs', 'Selected', fallback='')
            if selected_outputs_str:
                selected_names = selected_outputs_str.split(','); list_items = self.output_listbox.get(0, tk.END)
                self.output_listbox.selection_clear(0, tk.END)
                for name in selected_names:
                    if name in list_items:
                        try: idx = list_items.index(name); self.output_listbox.selection_set(idx)
                        except ValueError: pass
            self.show_isotherms_var.set(config.getboolean('PlotOptions', 'ShowIsotherms', fallback=False))
            self.show_isobars_var.set(config.getboolean('PlotOptions', 'ShowIsobars', fallback=False))
            if config.has_section('RankineInputs'):
                if config.has_option('RankineInputs', 'PBoiler'): self.rankine_p_boiler_var.set(config.get('RankineInputs', 'PBoiler'))
                if config.has_option('RankineInputs', 'TBoiler'): self.rankine_t_boiler_var.set(config.get('RankineInputs', 'TBoiler'))
                if config.has_option('RankineInputs', 'PCondenser'): self.rankine_p_cond_var.set(config.get('RankineInputs', 'PCondenser'))
                if config.has_option('RankineInputs', 'EtaTurbine'): self.rankine_eta_turbine_var.set(config.get('RankineInputs', 'EtaTurbine'))
        except Exception as e: logging.error(f"Error loading settings: {e}", exc_info=True); messagebox.showwarning("Load Settings Error", f"Could not load settings.\nUsing defaults.\n\nError: {e}", parent=self); self._setup_variables(); self._initial_state_setup()

    def _save_settings(self):
        logging.info(f"Saving settings to {CONFIG_FILE}..."); config = configparser.ConfigParser()
        try:
            config['Window'] = {}; config['Window']['Geometry'] = self.geometry()
            config['State'] = {}; config['State']['LastFluid'] = self.fluid_var.get()
            config['Inputs'] = {}; config['Inputs']['Prop1Name'] = self.prop1_name_var.get(); config['Inputs']['Prop1Value'] = self.prop1_value_var.get(); config['Inputs']['Prop1Unit'] = self.prop1_unit_var.get(); config['Inputs']['Prop2Name'] = self.prop2_name_var.get(); config['Inputs']['Prop2Value'] = self.prop2_value_var.get(); config['Inputs']['Prop2Unit'] = self.prop2_unit_var.get()
            config['Outputs'] = {}; selected_indices = self.output_listbox.curselection(); selected_names = [self.output_listbox.get(i) for i in selected_indices]; config['Outputs']['Selected'] = ','.join(selected_names)
            config['PlotOptions'] = {}; config['PlotOptions']['ShowIsotherms'] = str(self.show_isotherms_var.get()); config['PlotOptions']['ShowIsobars'] = str(self.show_isobars_var.get())
            config['RankineInputs'] = {}; config['RankineInputs']['PBoiler'] = self.rankine_p_boiler_var.get(); config['RankineInputs']['TBoiler'] = self.rankine_t_boiler_var.get(); config['RankineInputs']['PCondenser'] = self.rankine_p_cond_var.get(); config['RankineInputs']['EtaTurbine'] = self.rankine_eta_turbine_var.get()
            with open(CONFIG_FILE, 'w') as configfile: config.write(configfile)
            logging.info("Settings saved successfully.")
        except Exception as e: logging.error(f"Error saving settings: {e}", exc_info=True); messagebox.showerror("Save Settings Error", f"Could not save settings:\n{e}", parent=self)

    def _on_closing(self):
        logging.info("Closing application and saving settings..."); self._save_settings(); self.destroy()

    def _filter_fluid_list(self, event=None):
        filter_text = self.fluid_filter_var.get().lower(); self.fluid_listbox.delete(0, tk.END)
        for fluid in self._full_fluid_list:
            if filter_text in fluid.lower(): self.fluid_listbox.insert(tk.END, fluid)

    def _on_fluid_select_from_list(self, event=None):
        selected_indices = self.fluid_listbox.curselection()
        if not selected_indices: return
        selected_fluid = self.fluid_listbox.get(selected_indices[0])
        if selected_fluid:
            logging.info(f"Fluid selected via listbox: {selected_fluid}"); self.fluid_var.set(selected_fluid)
            if self.current_fluid_for_plot is None or self.current_fluid_for_plot != selected_fluid:
                self._update_fluid_info(); logging.info(f"Fluid changed. Clearing plots.")
                self._clear_ph_plot(); self._clear_ts_plot(); self.current_fluid_for_plot = None
                self._clear_rankine_results(); self.last_si_results = None; self._reset_output_display() ;self._clear_refrig_results()
                self._trigger_ph_plot_update(); self._trigger_ts_plot_update()
                self._update_psat_plot(selected_fluid) # NEW: Update P-T plot immediately
            else: self._update_fluid_info()

    def _get_prop_code_from_display(self, d): return self.display_to_code_map.get(d, 'DEFAULT')
    def _get_prop_details(self, c): return PROPERTY_INFO.get(c.upper(), PROPERTY_INFO['DEFAULT'])
    def _get_pint_unit_options(self, c): return PINT_UNIT_OPTIONS.get(self._get_prop_details(c).prop_type, PINT_UNIT_OPTIONS['Default'])

    def _update_unit_options(self, event=None):
        widget = event.widget if event else None; current_prop1_unit = self.prop1_unit_var.get(); current_prop2_unit = self.prop2_unit_var.get()
        p1_code = self._get_prop_code_from_display(self.prop1_name_var.get()); p1_list = self._get_pint_unit_options(p1_code); self.prop1_unit_combo['values'] = p1_list
        if widget == self.prop1_name_combo or not current_prop1_unit or current_prop1_unit not in p1_list:
            pref_def = 'degC' if p1_code == 'T' else ('bar' if p1_code == 'P' else ('dimensionless' if p1_code == 'Q' else (p1_list[0] if p1_list else '')))
            unit_to_set_p1 = pref_def if pref_def in p1_list else (p1_list[0] if p1_list else '')
        else: unit_to_set_p1 = current_prop1_unit
        self.prop1_unit_var.set(unit_to_set_p1)
        try: self.prop1_unit_combo.current(p1_list.index(unit_to_set_p1) if unit_to_set_p1 in p1_list else 0)
        except: pass

        p2_code = self._get_prop_code_from_display(self.prop2_name_var.get()); p2_list = self._get_pint_unit_options(p2_code); self.prop2_unit_combo['values'] = p2_list
        if widget == self.prop2_name_combo or not current_prop2_unit or current_prop2_unit not in p2_list:
            pref_def = 'bar' if p2_code == 'P' else ('degC' if p2_code == 'T' else ('dimensionless' if p2_code == 'Q' else (p2_list[0] if p2_list else '')))
            unit_to_set_p2 = pref_def if pref_def in p2_list else (p2_list[0] if p2_list else '')
        else: unit_to_set_p2 = current_prop2_unit
        self.prop2_unit_var.set(unit_to_set_p2)
        try: self.prop2_unit_combo.current(p2_list.index(unit_to_set_p2) if unit_to_set_p2 in p2_list else 0)
        except: pass

    def _update_fluid_info(self, event=None):
        fluid = self.fluid_var.get(); info_widgets = {"Tcrit": self.fluid_info_tcrit_var, "Pcrit": self.fluid_info_pcrit_var, "M": self.fluid_info_molmass_var, "Tmin": self.fluid_info_tmin_var, "Tmax": self.fluid_info_tmax_var}
        if not fluid:
            for tk_var in info_widgets.values(): tk_var.set(NA_TEXT)
            return
        fluid_info_si = self.calculator.get_fluid_info(fluid)
        for key, tk_var in info_widgets.items():
            val_si = fluid_info_si.get(key)
            if val_si == NA_TEXT or not isinstance(val_si, (int, float)): tk_var.set(NA_TEXT); continue
            try:
                prop_code = {'Tcrit':'T', 'Pcrit':'P', 'M':'M', 'Tmin':'T', 'Tmax':'T'}.get(key); si_unit_str = PROPERTY_INFO[prop_code].si_unit_str
                qty = self.Q_(val_si, si_unit_str); tk_var.set(f"{qty.to_compact():~P.6g}")
            except: tk_var.set(f"{val_si:.6g} (?)")

    def _clear_values(self):
        logging.info("Clearing input fields and results."); self.fluid_var.set("Water")
        self.prop1_name_var.set(self.code_to_display_map.get("T", "Temperature")); self.prop1_value_var.set("25")
        self.prop2_name_var.set(self.code_to_display_map.get("P", "Pressure")); self.prop2_value_var.set("1.01325")
        self._update_unit_options(); self.fluid_filter_var.set("Water"); self._filter_fluid_list()
        try:
            list_items = self.fluid_listbox.get(0, tk.END)
            if "Water" in list_items: idx = list_items.index("Water"); self.fluid_listbox.selection_clear(0, tk.END); self.fluid_listbox.selection_set(idx); self.fluid_listbox.see(idx)
        except: pass
        self._reset_output_display(); self.show_isotherms_var.set(False); self.show_isobars_var.set(False)
        self._clear_rankine_results();self._clear_refrig_results()
        if hasattr(self, 'psychro_tree'): [self.psychro_tree.delete(i) for i in self.psychro_tree.get_children()]
        self._clear_ph_plot(); self._clear_ts_plot(); self._clear_psat_plot()
        self.last_si_results = None; self.current_fluid_for_plot = None; self._update_fluid_info()
        self._update_psat_plot(self.fluid_var.get()) # NEW: Redraw P-T plot for Water
        if hasattr(self, 'copy_button'): self.copy_button.config(text=COPY_BUTTON_TEXT, state=tk.NORMAL)

    def _copy_result_to_clipboard(self):
        try:
            data = ["\t".join([self.result_tree.heading(c)['text'] for c in self.result_tree['columns']])]; general_error_item = self.result_tree.tag_has('generalerrorrow')
            for item_id in self.result_tree.get_children():
                item = self.result_tree.item(item_id); values = item['values']; tags = item['tags']
                if 'generalerrorrow' in tags: continue
                data.append("\t".join(map(str, values)))
            if general_error_item: error_values = self.result_tree.item(general_error_item[0])['values']; data.append("\n" + "\t".join(map(str, error_values)))
            clip_text = "\n".join(data)
            if len(data) <= 1 and not general_error_item : messagebox.showinfo("Copy Info", "No results to copy.", parent=self); return
            self.clipboard_clear(); self.clipboard_append(clip_text); logging.info(f"Copied {len(data)-1} result rows to clipboard.")
            self.copy_button.config(text=COPIED_BUTTON_TEXT, state=tk.DISABLED); self.copy_button.after(1500, lambda: self.copy_button.config(text=COPY_BUTTON_TEXT, state=tk.NORMAL))
        except Exception as e: logging.error(f"Copy to clipboard failed: {e}", exc_info=True); messagebox.showerror("Error", f"Could not copy results:\n{e}", parent=self)

    def _calculate_property_orchestrator(self):
        self._reset_output_display(); validated_si = self._validate_inputs()
        if validated_si is None: return
        logging.info("Input validation successful. Proceeding with calculation..."); indices = self.output_listbox.curselection()
        if not indices: messagebox.showinfo("Info", "Please select one or more output properties.", parent=self); return
        req_codes = [self._get_prop_code_from_display(self.output_listbox.get(i)) for i in indices]; req_codes = [c for c in req_codes if c != 'DEFAULT']
        if not req_codes: messagebox.showerror("Error", "Could not map selected outputs.", parent=self); return
        calc_codes = set(req_codes) | {'P', 'H', 'T', 'S'}; logging.info(f"Calculating properties for {validated_si['fluid']}...")
        try: self.last_si_results = self.calculator.calculate_properties(validated_si["fluid"], validated_si["prop1_code"], validated_si["prop1_value_si"], validated_si["prop2_code"], validated_si["prop2_value_si"], list(calc_codes))
        except Exception as e:
            logging.error(f"Unexpected error: {e}", exc_info=True); self._reset_output_display()
            self.result_tree.insert('', tk.END, values=("Fatal Error", f"Calculation failed: {e}", "-"), tags=('generalerrorrow',))
            self.last_si_results = None; self._clear_ph_plot(); self._clear_ts_plot(); self.current_fluid_for_plot = None; return
        self._update_output_display(self.last_si_results, req_codes); current_fluid = validated_si["fluid"]
        self._clear_cycle_lines_from_plots()
        if self.last_si_results and not self.last_si_results.get("error"):
            redraw_dome_needed = (current_fluid != self.current_fluid_for_plot)
            self._update_ph_plot(self.last_si_results, current_fluid, redraw_dome=redraw_dome_needed)
            self._update_ts_plot(self.last_si_results, current_fluid, redraw_dome=redraw_dome_needed)
            if redraw_dome_needed: self.current_fluid_for_plot = current_fluid
        elif self.last_si_results and self.last_si_results.get("error"):
             self._clear_ph_plot(); self._clear_ts_plot(); self.current_fluid_for_plot = None
        else: self._clear_ph_plot(); self._clear_ts_plot(); self.current_fluid_for_plot = None

    def _validate_inputs(self):
        fluid = self.fluid_var.get()
        if not fluid: return self._handle_validation_error("Please select a fluid.", "No fluid selected")
        p1_name = self.prop1_name_var.get(); p2_name = self.prop2_name_var.get(); p1_code = self._get_prop_code_from_display(p1_name); p2_code = self._get_prop_code_from_display(p2_name)
        p1_val_str = self.prop1_value_var.get().strip(); p2_val_str = self.prop2_value_var.get().strip(); p1_unit_str = self.prop1_unit_var.get(); p2_unit_str = self.prop2_unit_var.get()
        if p1_code=='DEFAULT' or p2_code=='DEFAULT': return self._handle_validation_error("Select valid input properties.", "Input prop missing/invalid")
        if not p1_val_str or not p2_val_str: return self._handle_validation_error("Enter values for both properties.", "Input val missing")
        if not p1_unit_str or not p2_unit_str: return self._handle_validation_error("Select units for both properties.", "Input unit missing")
        try: p1_val = float(p1_val_str); p2_val = float(p2_val_str)
        except ValueError: return self._handle_validation_error(f"Input values must be numeric.", "Non-numeric input")
        if p1_code == p2_code: return self._handle_validation_error("Input properties must be different.", "Identical inputs")
        try:
            q1 = self.Q_(p1_val, p1_unit_str); q2 = self.Q_(p2_val, p2_unit_str); p1_si_unit = self._get_prop_details(p1_code).si_unit_str; p2_si_unit = self._get_prop_details(p2_code).si_unit_str
            p1_si_mag = q1.to_base_units().magnitude if p1_code == 'Q' and p1_si_unit == 'dimensionless' else q1.to(p1_si_unit).magnitude
            p2_si_mag = q2.to_base_units().magnitude if p2_code == 'Q' and p2_si_unit == 'dimensionless' else q2.to(p2_si_unit).magnitude
        except Exception as e: return self._handle_validation_error(f"Unit conversion error: {e}", f"Pint Error: {e}")
        if p1_code == 'Q' or p2_code == 'Q':
            other_code = p2_code if p1_code == 'Q' else p1_code
            if other_code not in ['T', 'P', 'H', 'S', 'U']: return self._handle_validation_error("Quality (Q) must be paired with T, P, H, S, or U.", f"Invalid Q pair ({other_code})")
            q_val_si = p1_si_mag if p1_code == 'Q' else p2_si_mag
            if not (-0.00001 <= q_val_si <= 1.00001): return self._handle_validation_error(f"Quality (Q) must be between 0 and 1.", "Q range")
            q_val_si_clamped = max(0.0, min(1.0, q_val_si))
            if p1_code == 'Q': p1_si_mag = q_val_si_clamped
            else: p2_si_mag = q_val_si_clamped
        return {"fluid": fluid, "prop1_code": p1_code, "prop1_value_si": p1_si_mag, "prop2_code": p2_code, "prop2_value_si": p2_si_mag}

    def _handle_validation_error(self, user_msg, log_msg):
        logging.error(f"Input Validation Error: {log_msg}"); messagebox.showerror("Input Error", user_msg, parent=self); self._reset_output_display(); self.last_si_results = None; return None

    def _reset_output_display(self):
        if hasattr(self, 'result_tree'):
            for item in self.result_tree.get_children(): self.result_tree.delete(item)
        self.phase_var.set(DEFAULT_PHASE_TEXT)
        if hasattr(self, 'copy_button'): self.copy_button.config(text=COPY_BUTTON_TEXT, state=tk.NORMAL)
        self.update_idletasks()

    def _update_output_display(self, results_dict, requested_codes):
        self._reset_output_display()
        if not results_dict: self.result_tree.insert('', tk.END, values=("Error", "No results", "-"), tags=('generalerrorrow',)); return
        self.phase_var.set(f"{PHASE_PREFIX}{results_dict.get('phase', 'Unknown')}")
        general_error = results_dict.get("error")
        if general_error: self.result_tree.insert('', 0, values=("Calculation Status", general_error, "-"), tags=('generalerrorrow',))
        calc_vals = results_dict.get("values", {})
        display_order = []
        for code in requested_codes:
             if code in calc_vals and code not in display_order: display_order.append(code)
        for i, code in enumerate(display_order):
            value = calc_vals.get(code); prop_details = self._get_prop_details(code); is_error = False
            if isinstance(value, (int, float)):
                try:
                    qty = self.Q_(value, prop_details.si_unit_str)
                    val_str = f"{qty.magnitude:.5g}" if abs(qty.magnitude) > 1e6 or (abs(qty.magnitude) < 1e-3 and abs(qty.magnitude) != 0) else f"{qty.magnitude:.6f}".rstrip('0').rstrip('.')
                    if val_str == "": val_str = "0"
                    unit_disp = f"{qty.units:~P}"
                except: val_str = f"{value:.6g}"; unit_disp = prop_details.si_unit_str
            else: is_error = True; val_str = str(value) if value is not None else NA_TEXT; unit_disp = "-"
            tags = ('oddrow', 'errorrow') if is_error and (i%2==1) else (('errorrow',) if is_error else ('oddrow',) if (i+(1 if general_error else 0))%2==1 else ('evenrow',))
            self.result_tree.insert('', tk.END, values=(prop_details.description, val_str, unit_disp), tags=tags)
        self.result_tree.tag_configure('oddrow', background='#F0F0F0'); self.result_tree.tag_configure('evenrow', background='white'); self.result_tree.tag_configure('errorrow', foreground='red'); self.result_tree.tag_configure('generalerrorrow', foreground='red', font=('TkDefaultFont', 9, 'bold'))
        if general_error and not calc_vals and hasattr(self, 'copy_button'): self.copy_button.config(state=tk.DISABLED)

    def _trigger_ph_plot_update(self):
        if self.current_fluid_for_plot: self._update_ph_plot(self.last_si_results if self.last_si_results else {}, self.current_fluid_for_plot, redraw_dome=False, force_redraw_lines=True)
        else: self._clear_ph_plot()

    def _trigger_ts_plot_update(self):
        if self.current_fluid_for_plot: self._update_ts_plot(self.last_si_results if self.last_si_results else {}, self.current_fluid_for_plot, redraw_dome=False, force_redraw_lines=True)
        else: self._clear_ts_plot()

    def _plot_ph_base(self, fluid_name, dome_data):
        try:
            self.ph_ax.clear(); self.ph_ax.set_xlabel("Specific Enthalpy (h) [kJ/kg]"); self.ph_ax.set_ylabel("Pressure (P) [kPa]"); self.ph_ax.set_yscale('log')
            if dome_data.get('error'): self.ph_ax.set_title(f"P-h - {fluid_name} (Error)")
            elif dome_data.get('p') is not None:
                self.ph_ax.plot(dome_data['h_liq'], dome_data['p'], 'b-', lw=1.5, zorder=5); self.ph_ax.plot(dome_data['h_vap'], dome_data['p'], 'r-', lw=1.5, zorder=5)
                if dome_data.get('crit_point'): self.ph_ax.plot(dome_data['crit_point'][0], dome_data['crit_point'][1], 'ko', ms=5, zorder=10)
                self.ph_ax.set_title(f"P-h - {fluid_name}")
                try: self.ph_ax.set_xlim(min(dome_data['h_liq'])*0.9, max(dome_data['h_vap'])*1.1); self.ph_ax.set_ylim(min(dome_data['p'])*0.8, max(dome_data['p'])*1.2)
                except: pass
            else: self.ph_ax.set_title(f"P-h - {fluid_name} (No Data)")
            
            if self.show_isotherms_var.get() and self.last_si_results:
                t_si = self.last_si_results.get('values', {}).get('T')
                if isinstance(t_si, (int, float)):
                    iso = self.calculator.calculate_ph_isotherm(fluid_name, t_si)
                    if not iso.get('error'): self.ph_ax.plot(iso['h'], iso['p'], 'g--', lw=0.8, label=f'T={self.Q_(t_si,"K").to("degC").m:.0f}C', zorder=3); self.ph_ax.legend(loc='best', fontsize='small')
            self.ph_fig.tight_layout()
        except Exception: pass

    def _plot_ts_base(self, fluid_name, dome_data):
        try:
            self.ts_ax.clear(); self.ts_ax.set_xlabel("Entropy [kJ/(kg*K)]"); self.ts_ax.set_ylabel("Temperature [K]")
            if dome_data.get('error'): self.ts_ax.set_title(f"T-s - {fluid_name} (Error)")
            elif dome_data.get('t') is not None:
                self.ts_ax.plot(dome_data['s_liq'], dome_data['t'], 'b-', lw=1.5, zorder=5); self.ts_ax.plot(dome_data['s_vap'], dome_data['t'], 'r-', lw=1.5, zorder=5)
                if dome_data.get('crit_point'): self.ts_ax.plot(dome_data['crit_point'][0], dome_data['crit_point'][1], 'ko', ms=5, zorder=10)
                self.ts_ax.set_title(f"T-s - {fluid_name}")
                try: self.ts_ax.set_xlim(min(dome_data['s_liq'])*0.9, max(dome_data['s_vap'])*1.1); self.ts_ax.set_ylim(min(dome_data['t'])*0.9, max(dome_data['t'])*1.1)
                except: pass
            else: self.ts_ax.set_title(f"T-s - {fluid_name} (No Data)")

            if self.show_isobars_var.get() and self.last_si_results:
                p_si = self.last_si_results.get('values', {}).get('P')
                if isinstance(p_si, (int, float)):
                    iso = self.calculator.calculate_ts_isobar(fluid_name, p_si)
                    if not iso.get('error'): self.ts_ax.plot(iso['s'], iso['t'], 'm--', lw=0.8, label=f'P={self.Q_(p_si,"Pa").to("bar").m:.2f}bar', zorder=3); self.ts_ax.legend(loc='best', fontsize='small')
            self.ts_fig.tight_layout()
        except Exception: pass

    def _clear_ph_plot(self):
        try:
            if hasattr(self, 'ph_ax'): self.ph_ax.clear(); self.ph_ax.set_xlabel("Specific Enthalpy (h) [kJ/kg]"); self.ph_ax.set_ylabel("Pressure (P) [kPa]"); self.ph_ax.set_yscale('log'); self.ph_ax.set_title("P-h Diagram"); self.ph_ax.set_xlim(0, 1); self.ph_ax.set_ylim(1, 10); self.ph_fig.tight_layout(); self.ph_canvas.draw_idle()
        except: pass

    def _update_ph_plot(self, si_results, fluid_name, redraw_dome, force_redraw_lines=False):
        try:
            calc_vals = si_results.get("values", {}) if si_results else {}
            h_plot = self.Q_(calc_vals.get('H'),'J/kg').to('kJ/kg').m if isinstance(calc_vals.get('H'),(int,float)) else None
            p_plot = self.Q_(calc_vals.get('P'),'pascal').to('kPa').m if isinstance(calc_vals.get('P'),(int,float)) else None
            if redraw_dome or force_redraw_lines:
                self._plot_ph_base(fluid_name, self.calculator.calculate_ph_dome(fluid_name))
            else:
                for item in [i for i in (self.ph_ax.lines + self.ph_ax.collections + self.ph_ax.texts) if i.get_label() in ('state_point_marker', 'state_point_annot')]: item.remove()
            if h_plot is not None and p_plot is not None:
                self.ph_ax.plot(h_plot, p_plot, 'X', c='lime', ms=10, mec='k', label='state_point_marker', zorder=15)
                self.ph_ax.annotate(si_results.get('phase','?'), (h_plot, p_plot), xytext=(7,-7), textcoords='offset points', size='small', label='state_point_annot', bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.6, ec="none"), zorder=16)
                try:
                     curr_xlim = self.ph_ax.get_xlim(); curr_ylim = self.ph_ax.get_ylim()
                     if h_plot < curr_xlim[0] or h_plot > curr_xlim[1] or p_plot < curr_ylim[0] or p_plot > curr_ylim[1]:
                         self.ph_ax.set_xlim(min(curr_xlim[0], h_plot*0.9), max(curr_xlim[1], h_plot*1.1))
                         self.ph_ax.set_ylim(min(curr_ylim[0], p_plot*0.5), max(curr_ylim[1], p_plot*1.5)); self.ph_fig.tight_layout()
                except: pass
            self.ph_canvas.draw_idle()
        except: pass

    def _clear_ts_plot(self):
        try:
            if hasattr(self, 'ts_ax'): self.ts_ax.clear(); self.ts_ax.set_xlabel("Entropy [kJ/(kg*K)]"); self.ts_ax.set_ylabel("Temperature [K]"); self.ts_ax.set_title("T-s Diagram"); self.ts_ax.set_xlim(0, 1); self.ts_ax.set_ylim(200, 300); self.ts_fig.tight_layout(); self.ts_canvas.draw_idle()
        except: pass

    def _update_ts_plot(self, si_results, fluid_name, redraw_dome, force_redraw_lines=False):
        try:
            calc_vals = si_results.get("values", {}) if si_results else {}
            s_plot = self.Q_(calc_vals.get('S'),'J/(kg*K)').to('kJ/(kg*K)').m if isinstance(calc_vals.get('S'),(int,float)) else None
            t_plot = self.Q_(calc_vals.get('T'),'K').to('K').m if isinstance(calc_vals.get('T'),(int,float)) else None
            if redraw_dome or force_redraw_lines:
                self._plot_ts_base(fluid_name, self.calculator.calculate_ts_dome(fluid_name))
            else:
                for item in [i for i in (self.ts_ax.lines + self.ts_ax.collections + self.ts_ax.texts) if i.get_label() in ('state_point_marker_ts', 'state_point_annot_ts')]: item.remove()
            if s_plot is not None and t_plot is not None:
                self.ts_ax.plot(s_plot, t_plot, 'X', c='lime', ms=10, mec='k', label='state_point_marker_ts', zorder=15)
                self.ts_ax.annotate(si_results.get('phase','?'), (s_plot, t_plot), xytext=(7,-7), textcoords='offset points', size='small', label='state_point_annot_ts', bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.6, ec="none"), zorder=16)
                try:
                    curr_xlim = self.ts_ax.get_xlim(); curr_ylim = self.ts_ax.get_ylim()
                    if s_plot < curr_xlim[0] or s_plot > curr_xlim[1] or t_plot < curr_ylim[0] or t_plot > curr_ylim[1]:
                         self.ts_ax.set_xlim(min(curr_xlim[0], s_plot*0.9), max(curr_xlim[1], s_plot*1.1))
                         self.ts_ax.set_ylim(min(curr_ylim[0], t_plot*0.9), max(curr_ylim[1], t_plot*1.1)); self.ts_fig.tight_layout()
                except: pass
            self.ts_canvas.draw_idle()
        except: pass

    # --- NEW: P-T Saturation Plot Methods ---
    def _clear_psat_plot(self):
        try:
            if hasattr(self, 'psat_ax'):
                self.psat_ax.clear(); self.psat_ax.set_xlabel("Temperature (T) [K]"); self.psat_ax.set_ylabel("Saturation Pressure (P) [kPa]"); self.psat_ax.set_yscale('log'); self.psat_ax.set_title("P-T Saturation Curve"); self.psat_ax.set_xlim(100, 200); self.psat_ax.set_ylim(1, 10); self.psat_fig.tight_layout(); self.psat_canvas.draw_idle()
        except: pass

    def _update_psat_plot(self, fluid_name):
        if not fluid_name or not hasattr(self, 'psat_ax'): self._clear_psat_plot(); return
        logging.info(f"Updating P-T Saturation plot for {fluid_name}")
        self._clear_psat_plot()
        try:
            psat_data = self.calculator.calculate_psat_vs_t(fluid_name)
            if psat_data.get('error'): self.psat_ax.set_title(f"P-T Saturation - {fluid_name} (Error)")
            elif psat_data.get('t') is not None:
                self.psat_ax.plot(psat_data['t'], psat_data['p'], 'b-', lw=1.5, zorder=5)
                self.psat_ax.set_title(f"P-T Saturation - {fluid_name}")
                try:
                    t_min, t_max = min(psat_data['t']), max(psat_data['t'])
                    p_min, p_max = min(psat_data['p']), max(psat_data['p'])
                    self.psat_ax.set_xlim(t_min * 0.95, t_max * 1.05)
                    self.psat_ax.set_ylim(p_min * 0.9, p_max * 1.1)
                except: pass
            else: self.psat_ax.set_title(f"P-T Saturation - {fluid_name} (No Data)")
            self.psat_fig.tight_layout(); self.psat_canvas.draw_idle()
        except Exception as e: logging.error(f"Error updating P-T plot: {e}")

    def _clear_rankine_results(self):
        if hasattr(self, 'rankine_tree'):
            for item in self.rankine_tree.get_children(): self.rankine_tree.delete(item)
        if hasattr(self, 'rankine_wpump_label'):
             for lbl in [self.rankine_wpump_label, self.rankine_qboiler_label, self.rankine_wturbine_label, self.rankine_qcond_label]: lbl.config(text=f"{lbl.cget('text').split(':')[0]}: -")
             self.rankine_eta_label.config(text="η_thermal: -")
        self._clear_cycle_lines_from_plots()

    def _clear_cycle_lines_from_plots(self):
        for ax in [self.ph_ax, self.ts_ax]:
             for item in [i for i in (ax.lines + ax.collections + ax.texts) if i.get_label() in ('cycle_point', 'cycle_line', 'cycle_label')]: item.remove()
        try: self.ph_canvas.draw_idle(); self.ts_canvas.draw_idle()
        except: pass

    def _calculate_rankine_orchestrator(self):
        self._clear_rankine_results(); fluid = self.fluid_var.get()
        if not fluid: messagebox.showerror("Input Error", "Select a fluid first."); return
        try:
            p_boiler = self.Q_(float(self.rankine_p_boiler_var.get()), 'kPa').to('Pa').m
            t_boiler = self.Q_(float(self.rankine_t_boiler_var.get()), 'degC').to('K').m
            p_cond = self.Q_(float(self.rankine_p_cond_var.get()), 'kPa').to('Pa').m
            eta = float(self.rankine_eta_turbine_var.get()) / 100.0
            if fluid != self.current_fluid_for_plot: self.current_fluid_for_plot = fluid
            cycle_results = self.calculator.calculate_rankine_cycle(fluid, p_boiler, t_boiler, p_cond, eta)
            if cycle_results.get('error'): messagebox.showerror("Calculation Error", cycle_results['error']); self.rankine_eta_label.config(text=f"Error: {cycle_results['error']}"); return
            states = cycle_results.get('states', {})
            for i in sorted(states.keys()):
                p = states[i]; q = p.get('Q', np.nan)
                vals = (i, f"{p.get('T',np.nan):.2f}", f"{self.Q_(p.get('P',np.nan),'Pa').to('kPa').m:.1f}", f"{self.Q_(p.get('H',np.nan),'J/kg').to('kJ/kg').m:.2f}", f"{self.Q_(p.get('S',np.nan),'J/(kg*K)').to('kJ/(kg*K)').m:.4f}", f"{q:.4f}" if 0<=q<=1 else ("-" if isnan(q) else (">1" if q>1 else "<0")))
                self.rankine_tree.insert('', tk.END, values=vals)
            m = cycle_results.get('metrics', {})
            self.rankine_wpump_label.config(text=f"W_pump: {self.Q_(m.get('W_pump',0),'J/kg').to('kJ/kg').m:.2f} kJ/kg")
            self.rankine_qboiler_label.config(text=f"Q_boiler: {self.Q_(m.get('Q_boiler',0),'J/kg').to('kJ/kg').m:.2f} kJ/kg")
            self.rankine_wturbine_label.config(text=f"W_turbine: {self.Q_(m.get('W_turbine',0),'J/kg').to('kJ/kg').m:.2f} kJ/kg")
            self.rankine_qcond_label.config(text=f"Q_condenser: {self.Q_(m.get('Q_condenser',0),'J/kg').to('kJ/kg').m:.2f} kJ/kg")
            self.rankine_eta_label.config(text=f"η_thermal: {m.get('eta_thermal',0)*100:.2f} %")
            self._update_cycle_plots(cycle_results, fluid, redraw_base=True)
        except Exception as e: messagebox.showerror("Error", f"Rankine error:\n{e}"); logging.error(f"Rankine error: {e}")

    def _update_cycle_plots(self, cycle_results, fluid_name, redraw_base):
        states = cycle_results.get('states');
        if not states: return
        try:
            pts = sorted(states.keys()); h_p, p_p, s_p, t_p = [], [], [], []
            for i in pts:
                s = states[i]; h_p.append(self.Q_(s.get('H',np.nan),'J/kg').to('kJ/kg').m); p_p.append(self.Q_(s.get('P',np.nan),'Pa').to('kPa').m); t_p.append(s.get('T',np.nan)); s_p.append(self.Q_(s.get('S',np.nan),'J/(kg*K)').to('kJ/(kg*K)').m)
            if any(isnan(v) for v in h_p+p_p+s_p+t_p): return
            self._clear_cycle_lines_from_plots()
            if redraw_base: self._plot_ph_base(fluid_name, self.calculator.calculate_ph_dome(fluid_name)); self._plot_ts_base(fluid_name, self.calculator.calculate_ts_dome(fluid_name))
            self.ph_ax.plot(h_p+[h_p[0]], p_p+[p_p[0]], 'r-', lw=1.0, label='cycle_line', zorder=11); self.ts_ax.plot(s_p+[s_p[0]], t_p+[t_p[0]], 'r-', lw=1.0, label='cycle_line', zorder=11)
            for i, pt in enumerate(pts):
                self.ph_ax.plot(h_p[i], p_p[i], 'ro', ms=5, label='cycle_point', zorder=12); self.ph_ax.text(h_p[i], p_p[i], f' {pt}', color='red', fontsize=9, label='cycle_label', zorder=13)
                self.ts_ax.plot(s_p[i], t_p[i], 'ro', ms=5, label='cycle_point', zorder=12); self.ts_ax.text(s_p[i], t_p[i], f' {pt}', color='red', fontsize=9, label='cycle_label', zorder=13)
            self.ph_canvas.draw_idle(); self.ts_canvas.draw_idle()
        except Exception as e: logging.error(f"Cycle plot error: {e}")
    # --- NEW: Refrigeration Cycle UI Methods ---

    def _create_refrig_widgets(self, parent_frame):
        """ Creates widgets for the Refrigeration Cycle tab. """
        # Input Frame
        input_frame = ttk.LabelFrame(parent_frame, text="Refrigeration Cycle Inputs", padding="10")
        input_frame.grid(row=0, column=0, padx=5, pady=5, sticky="new")

        # Row 0: Evap Pressure
        ttk.Label(input_frame, text="Evaporator Pressure (P_low):").grid(row=0, column=0, padx=5, pady=3, sticky="w")
        ttk.Entry(input_frame, textvariable=self.refrig_p_evap_var, width=10).grid(row=0, column=1, padx=5, pady=3, sticky="ew")
        ttk.Label(input_frame, text="kPa").grid(row=0, column=2, padx=(0,5), pady=3, sticky="w")

        # Row 1: Cond Pressure
        ttk.Label(input_frame, text="Condenser Pressure (P_high):").grid(row=1, column=0, padx=5, pady=3, sticky="w")
        ttk.Entry(input_frame, textvariable=self.refrig_p_cond_var, width=10).grid(row=1, column=1, padx=5, pady=3, sticky="ew")
        ttk.Label(input_frame, text="kPa").grid(row=1, column=2, padx=(0,5), pady=3, sticky="w")

        # Row 2: Superheat
        ttk.Label(input_frame, text="Superheat (ΔT):").grid(row=2, column=0, padx=5, pady=3, sticky="w")
        ttk.Entry(input_frame, textvariable=self.refrig_sh_var, width=10).grid(row=2, column=1, padx=5, pady=3, sticky="ew")
        ttk.Label(input_frame, text="K").grid(row=2, column=2, padx=(0,5), pady=3, sticky="w")

        # Row 3: Subcooling
        ttk.Label(input_frame, text="Subcooling (ΔT):").grid(row=3, column=0, padx=5, pady=3, sticky="w")
        ttk.Entry(input_frame, textvariable=self.refrig_sc_var, width=10).grid(row=3, column=1, padx=5, pady=3, sticky="ew")
        ttk.Label(input_frame, text="K").grid(row=3, column=2, padx=(0,5), pady=3, sticky="w")

        # Row 4: Efficiency
        ttk.Label(input_frame, text="Compressor Efficiency (η):").grid(row=4, column=0, padx=5, pady=3, sticky="w")
        ttk.Entry(input_frame, textvariable=self.refrig_eta_var, width=10).grid(row=4, column=1, padx=5, pady=3, sticky="ew")
        ttk.Label(input_frame, text="%").grid(row=4, column=2, padx=(0,5), pady=3, sticky="w")

        calc_button = ttk.Button(input_frame, text="Calculate Refrig Cycle", command=self._calculate_refrig_orchestrator)
        calc_button.grid(row=5, column=0, columnspan=3, pady=10)

        # Results Frame
        results_frame = ttk.LabelFrame(parent_frame, text="Cycle Performance", padding="10")
        results_frame.grid(row=0, column=1, rowspan=2, padx=5, pady=5, sticky="nsew")
        results_frame.grid_rowconfigure(1, weight=1); results_frame.grid_columnconfigure(0, weight=1)

        metrics_frame = ttk.Frame(results_frame, padding=(0, 5, 0, 10))
        metrics_frame.grid(row=0, column=0, sticky="new")

        self.refrig_copc_label = ttk.Label(metrics_frame, text="COP (Cooling): -", font="-weight bold")
        self.refrig_copc_label.grid(row=0, column=0, padx=5, pady=2, sticky="w")
        self.refrig_coph_label = ttk.Label(metrics_frame, text="COP (Heating): -")
        self.refrig_coph_label.grid(row=0, column=1, padx=15, pady=2, sticky="w")
        
        self.refrig_qevap_label = ttk.Label(metrics_frame, text="Q_evap: -")
        self.refrig_qevap_label.grid(row=1, column=0, padx=5, pady=1, sticky="w")
        self.refrig_wcomp_label = ttk.Label(metrics_frame, text="W_comp: -")
        self.refrig_wcomp_label.grid(row=1, column=1, padx=15, pady=1, sticky="w")

        # Treeview
        cols = ('point', 'T', 'P', 'H', 'S', 'Q')
        self.refrig_tree = ttk.Treeview(results_frame, columns=cols, show='headings', height=5)
        self.refrig_tree.grid(row=1, column=0, sticky="nsew")
        
        self.refrig_tree.heading('point', text='State')
        self.refrig_tree.heading('T', text='T [°C]') # Display Celsius for readability in HVAC
        self.refrig_tree.heading('P', text='P [kPa]')
        self.refrig_tree.heading('H', text='H [kJ/kg]')
        self.refrig_tree.heading('S', text='S [kJ/kg·K]')
        self.refrig_tree.heading('Q', text='Quality')
        
        for col in cols:
            width = 50 if col == 'point' else 85
            self.refrig_tree.column(col, width=width, anchor=tk.E if col != 'point' else tk.CENTER)
            
        scroll = ttk.Scrollbar(results_frame, orient=tk.VERTICAL, command=self.refrig_tree.yview)
        self.refrig_tree.configure(yscrollcommand=scroll.set)
        scroll.grid(row=1, column=1, sticky='ns')

    def _clear_refrig_results(self):
        if hasattr(self, 'refrig_tree'):
            for item in self.refrig_tree.get_children(): self.refrig_tree.delete(item)
        if hasattr(self, 'refrig_copc_label'):
            self.refrig_copc_label.config(text="COP (Cooling): -")
            self.refrig_coph_label.config(text="COP (Heating): -")
            self.refrig_qevap_label.config(text="Q_evap: -")
            self.refrig_wcomp_label.config(text="W_comp: -")
        # Note: Plots are cleared by the generic _clear_cycle_lines_from_plots called in orchestrator

    def _calculate_refrig_orchestrator(self):
        self._clear_refrig_results()
        fluid = self.fluid_var.get()
        if not fluid: messagebox.showerror("Input Error", "Select a fluid on the Calculator tab."); return

        try:
            # Get Inputs
            p_evap_kpa = float(self.refrig_p_evap_var.get())
            p_cond_kpa = float(self.refrig_p_cond_var.get())
            sh_k = float(self.refrig_sh_var.get())
            sc_k = float(self.refrig_sc_var.get())
            eta = float(self.refrig_eta_var.get()) / 100.0

            # Convert
            p_evap_pa = self.Q_(p_evap_kpa, 'kPa').to('Pa').m
            p_cond_pa = self.Q_(p_cond_kpa, 'kPa').to('Pa').m

            # Ensure plot consistency
            if fluid != self.current_fluid_for_plot: self.current_fluid_for_plot = fluid

            # Calc
            res = self.calculator.calculate_refrigeration_cycle(fluid, p_evap_pa, p_cond_pa, sh_k, sc_k, eta)

            if res.get('error'):
                messagebox.showerror("Error", res['error'])
                self.refrig_copc_label.config(text="Error")
                return

            # Display Metrics
            m = res.get('metrics', {})
            self.refrig_copc_label.config(text=f"COP (Cooling): {m.get('COP_cooling',0):.2f}")
            self.refrig_coph_label.config(text=f"COP (Heating): {m.get('COP_heating',0):.2f}")
            q_evap_kj = self.Q_(m.get('Q_evap',0), 'J/kg').to('kJ/kg').m
            w_comp_kj = self.Q_(m.get('W_comp',0), 'J/kg').to('kJ/kg').m
            self.refrig_qevap_label.config(text=f"Q_evap: {q_evap_kj:.1f} kJ/kg")
            self.refrig_wcomp_label.config(text=f"W_comp: {w_comp_kj:.1f} kJ/kg")

            # Display States
            states = res.get('states', {})
            for i in sorted(states.keys()):
                s = states[i]
                t_c = self.Q_(s['T'], 'K').to('degC').m
                p_kpa = self.Q_(s['P'], 'Pa').to('kPa').m
                h_kj = self.Q_(s['H'], 'J/kg').to('kJ/kg').m
                s_kj = self.Q_(s['S'], 'J/(kg*K)').to('kJ/(kg*K)').m
                q_val = s.get('Q', -1)
                q_str = f"{q_val:.2f}" if 0<=q_val<=1 else ("Sup" if q_val>1 else "Sub")
                
                self.refrig_tree.insert('', tk.END, values=(i, f"{t_c:.2f}", f"{p_kpa:.1f}", f"{h_kj:.1f}", f"{s_kj:.4f}", q_str))

            # Plot
            self._update_cycle_plots(res, fluid, redraw_base=True)
            
            # Switch to P-h tab automatically for better UX
            self.notebook.select(self.ph_plot_tab)

        except ValueError as ve: messagebox.showerror("Input Error", f"Invalid Input:\n{ve}")
        except Exception as e: messagebox.showerror("Error", f"Calculation failed:\n{e}"); logging.error(e, exc_info=True)
    # --- NEW: Psychrometric Tab Methods ---

    def _create_psychro_widgets(self, parent_frame):
        """ Creates widgets for the Psychrometric Calculator tab. """
        # --- Input Frame ---
        input_frame = ttk.LabelFrame(parent_frame, text="Air State Inputs", padding="10")
        input_frame.grid(row=0, column=0, padx=5, pady=5, sticky="new")
        input_frame.grid_columnconfigure(1, weight=1)

        # Dry Bulb Temp
        ttk.Label(input_frame, text="Dry Bulb Temp (Tdb):").grid(row=0, column=0, padx=5, pady=3, sticky="w")
        ttk.Entry(input_frame, textvariable=self.psychro_tdb_var, width=10).grid(row=0, column=1, padx=5, pady=3, sticky="ew")
        ttk.Label(input_frame, text="°C").grid(row=0, column=2, padx=(0,5), pady=3, sticky="w")

        # Relative Humidity
        ttk.Label(input_frame, text="Relative Humidity (RH):").grid(row=1, column=0, padx=5, pady=3, sticky="w")
        ttk.Entry(input_frame, textvariable=self.psychro_rh_var, width=10).grid(row=1, column=1, padx=5, pady=3, sticky="ew")
        ttk.Label(input_frame, text="%").grid(row=1, column=2, padx=(0,5), pady=3, sticky="w")

        # Atmospheric Pressure
        ttk.Label(input_frame, text="Atm. Pressure (Patm):").grid(row=2, column=0, padx=5, pady=3, sticky="w")
        ttk.Entry(input_frame, textvariable=self.psychro_p_var, width=10).grid(row=2, column=1, padx=5, pady=3, sticky="ew")
        ttk.Label(input_frame, text="kPa").grid(row=2, column=2, padx=(0,5), pady=3, sticky="w")

        # Calculate Button
        calc_btn = ttk.Button(input_frame, text="Calculate Properties", command=self._calculate_psychro_orchestrator)
        calc_btn.grid(row=3, column=0, columnspan=3, pady=10)

        # --- Results Frame ---
        results_frame = ttk.LabelFrame(parent_frame, text="Humid Air Properties", padding="10")
        results_frame.grid(row=0, column=1, padx=5, pady=5, sticky="nsew")
        results_frame.grid_rowconfigure(0, weight=1); results_frame.grid_columnconfigure(0, weight=1)

        cols = ('prop', 'val', 'unit')
        self.psychro_tree = ttk.Treeview(results_frame, columns=cols, show='headings', height=8)
        self.psychro_tree.grid(row=0, column=0, sticky="nsew")

        self.psychro_tree.heading('prop', text='Property')
        self.psychro_tree.heading('val', text='Value')
        self.psychro_tree.heading('unit', text='Unit')

        self.psychro_tree.column('prop', width=150, anchor=tk.W)
        self.psychro_tree.column('val', width=100, anchor=tk.E)
        self.psychro_tree.column('unit', width=100, anchor=tk.W)

        scroll = ttk.Scrollbar(results_frame, orient=tk.VERTICAL, command=self.psychro_tree.yview)
        self.psychro_tree.configure(yscrollcommand=scroll.set)
        scroll.grid(row=0, column=1, sticky='ns')

    def _calculate_psychro_orchestrator(self):
        # Clear previous
        for item in self.psychro_tree.get_children(): self.psychro_tree.delete(item)
        
        try:
            # Get & Convert Inputs
            t_db_c = float(self.psychro_tdb_var.get())
            rh_percent = float(self.psychro_rh_var.get())
            p_kpa = float(self.psychro_p_var.get())

            t_db_k = self.Q_(t_db_c, 'degC').to('K').m
            rh_frac = rh_percent / 100.0
            p_pa = self.Q_(p_kpa, 'kPa').to('Pa').m

            # Calculate
            result = self.calculator.calculate_humid_air_properties(t_db_k, rh_frac, p_pa)

            if result.get('error'):
                messagebox.showerror("Psychro Error", result['error'], parent=self)
                return

            # Display Results (Unit Conversions for readability)
            props = result.get('properties', {})
            
            # Define display logic: Key -> (Description, DisplayUnitPintStr, DisplayFormat)
            display_map = {
                'Tdp': ('Dew Point Temp', 'degC', '.2f'),
                'Twb': ('Wet Bulb Temp', 'degC', '.2f'),
                'W':   ('Humidity Ratio', 'dimensionless', '.5f'), # kg_w/kg_da
                'H':   ('Specific Enthalpy', 'kJ/kg', '.2f'),      # per kg dry air
                'V':   ('Specific Volume', 'm**3/kg', '.4f'),      # per kg dry air
                'M':   ('Viscosity', 'uPa*s', '.2f'),
                'K':   ('Conductivity', 'W/(m*K)', '.4f')
            }

            for key, val in props.items():
                if val is None: continue
                desc, unit_str, fmt = display_map.get(key, (key, '', '.4f'))
                
                # Conversion handling
                try:
                    # SI Unit coming from CoolProp
                    si_unit = 'K' if key in ['Tdp', 'Twb'] else ('J/kg' if key=='H' else ('m**3/kg' if key=='V' else ('Pa*s' if key=='M' else ('W/(m*K)' if key=='K' else 'dimensionless'))))
                    
                    qty = self.Q_(val, si_unit)
                    
                    if key == 'W':
                        # Special case for W: show as kg/kg (dimensionless) or g/kg? 
                        # Let's show kg/kg for standard, maybe append g/kg in text
                        disp_val = f"{val:{fmt}}"
                        disp_unit = "kg_w/kg_da"
                        # Add g/kg row ? optional.
                    else:
                        qty_conv = qty.to(unit_str)
                        disp_val = f"{qty_conv.magnitude:{fmt}}"
                        disp_unit = f"{qty_conv.units:~P}"
                    
                    self.psychro_tree.insert('', tk.END, values=(desc, disp_val, disp_unit))
                    
                except Exception as e:
                    logging.error(f"Formatting error for {key}: {e}")
                    self.psychro_tree.insert('', tk.END, values=(desc, str(val), "Error"))

        except ValueError as ve:
            messagebox.showerror("Input Error", f"Invalid numeric input: {ve}", parent=self)
        except Exception as e:
            messagebox.showerror("Error", f"Calculation failed: {e}", parent=self)
            logging.error(f"Psychro Orchestrator Error: {e}", exc_info=True)
if __name__ == "__main__":
    try: from ctypes import windll; windll.shcore.SetProcessDpiAwareness(1)
    except: pass
    try: app = CoolPropApp(); app.mainloop()
    except Exception as main_e:
        logging.critical(f"Fatal error: {main_e}", exc_info=True)
        try: root = tk.Tk(); root.withdraw(); messagebox.showerror("Fatal Application Error", f"Critical error:\n\n{main_e}"); root.destroy()
        except: pass
        sys.exit(1)