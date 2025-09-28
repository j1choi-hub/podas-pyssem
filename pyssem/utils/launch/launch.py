from sympy import zeros, Matrix, symbols
import pandas as pd
from datetime import datetime, timedelta
import numpy as np
from tqdm import tqdm
import os
import scipy.stats

def find_alt_bin(altitude, scen_properties):
    # Convert altitude ranges to numpy arrays for vectorized operations
    lower = np.array(scen_properties['R02'][:-1])
    upper = np.array(scen_properties['R02'][1:])
    
    # Create a boolean array where True indicates altitude is within the shell bounds
    shell_logic = (lower < altitude) & (altitude <= upper)
    
    # Find the index (or indices) where shell_logic is True
    shell_indices = np.where(shell_logic)[0]
    
    # Return the first index if found, otherwise return NaN
    shell_index = shell_indices[0] + 1 if shell_indices.size > 0 else np.nan
    return shell_index

def launch_func_null(t, h, species_properties, scen_properties):
    """
    No launch function for species without a launch function.

    Args:
        t (float): Time from scenario start in years
        h (array_like): The set of altitudes of the scenario above ellipsoid in km of shell lower edges.
        species_properties (dict): A dictionary with properties for the species
        scen_properties (dict): A dictionary with properties for the scenario

    Returns:
        numpy.ndarray: The rate of change in the species in each shell at the specified time due to launch.
                       If only one value is applied, it is assumed to be true for all shells.
    """
    if len(h) != scen_properties.n_shells:
        raise ValueError("Constant launch rate must be specified per altitude shell.")

    n_shells = scen_properties.n_shells
    scen_times = scen_properties.scen_times
    simulation_duration = scen_properties.simulation_duration
    # Calculate the number of time steps per year
    time_steps_per_year = len(scen_times) / simulation_duration

    # Create a matrix filled with zeros
    launch_rate_matrix = np.zeros((n_shells, len(scen_times)))

    return launch_rate_matrix

    # Lambdadot = zeros(scen_properties.n_shells, 1)
    #
    # for k in range(scen_properties.n_shells):
    #     Lambdadot[k, 0] = 0 * species_properties.sym[k]
    #
    # Lambdadot_list = [Lambdadot[k, 0] for k in range(scen_properties.n_shells)]
    #
    # return Lambdadot_list

def launch_func_constant(t, h, species_properties, scen_properties):
    """
    Adds a constant launch rate from species_properties.lambda_constant.
    Given a certain altitude, this will return the rate of change in the species in each shell at the specified time due to launch.

    Args:
        t (float): Time from scenario start in years.
        h (list or numpy.ndarray): Altitudes of the scenario above ellipsoid in km of shell lower edges.
        species_properties (dict): Properties for the species, including 'lambda_constant'.
        scen_properties (dict): Properties for the scenario, including 'N_shell'.

    Returns:
        list: np., a list of symbolic expressions representing the rate of change in the species in each shell due to launch.
    """
    if len(h) != scen_properties.n_shells:
        raise ValueError("Constant launch rate must be specified per altitude shell.")

    n_shells = scen_properties.n_shells
    scen_times = scen_properties.scen_times
    simulation_duration = scen_properties.simulation_duration
    # Calculate the number of time steps per year
    time_steps_per_year = len(scen_times) / simulation_duration

    # Create a matrix filled with zeros
    launch_rate_matrix = np.zeros((n_shells, len(scen_times)))

    # Find the index of the shell closest to the specified altitude
    for alt in scen_properties.R0_km:
        h_inds = np.argmin(np.abs(h - alt))

        # Calculate the number of satellites per time step
        launch_rate_per_year = species_properties.lambda_constant  # satellites/year
        launch_rate_per_time_step = launch_rate_per_year / time_steps_per_year

        # Fill the row corresponding to the nearest altitude with the launch rate per time step
        launch_rate_matrix[h_inds, :] = launch_rate_per_time_step

    return launch_rate_matrix


    # if len(h) != scen_properties.n_shells:
    #     raise ValueError("Constant launch rate must be specified per altitude shell.")
    #
    # # Create a symbolic variable for the launch rate
    # lambda_constant = symbols('lambda_constant')
    #
    # # Assign the constant launch rate to each shell
    # Lambdadot = Matrix(scen_properties.n_shells, 1, lambda i, j: lambda_constant)
    #
    # # Convert the Matrix of symbolic expressions to a list
    # Lambdadot_list = [Lambdadot[i] for i in range(scen_properties.n_shells)]
    #
    # print(Lambdadot_list)
    #
    # return Lambdadot_list

def launch_func_lambda_fun(t, h, species_properties, scen_properties):
    """
    This function will return the lambda function for a required species. 

    :param t: The time from the scenario start in years
    :type t: int
    :param h: The altitude above the ellipsoid in km of shell lower edge
    :type h: int
    :param species_properties: Species properties
    :type species_properties: Species
    :param scen_properties: Scenario properties
    :type scen_properties: ScenarioProperties
    :return: Lambdadot is the rate of change in the species in each sheel at the specified time due to launch
    :rtype: SciPy interp1d function
    """
    # # Find the index for the given altitude
    # h_inds = np.where(scen_properties.HMid == h)
    # print(species_properties.sym_name)

    # # Retrieve the appropriate lambda function for the altitude and evaluate it at time t
    # Lambdadot = species_properties.lambda_funs
    # return Lambdadot
    if len(h) != scen_properties.n_shells:
        raise ValueError("Constant launch rate must be specified per altitude shell.")

    mu_lam = 510
    sig_lam = 150**2
    #6000
    lam0 = 4000*(scipy.stats.norm.cdf(max(scen_properties.R0_km), loc=mu_lam, scale=np.sqrt(sig_lam)) - scipy.stats.norm.cdf(min(scen_properties.R0_km), loc=mu_lam, scale=np.sqrt(sig_lam)))

    n_shells = scen_properties.n_shells
    scen_times = scen_properties.scen_times
    simulation_duration = scen_properties.simulation_duration
    # Calculate the number of time steps per year
    time_steps_per_year = scen_properties.steps/simulation_duration

    # Create a matrix filled with zeros
    launch_rate_matrix = np.zeros((n_shells, len(scen_times)))
    #print("R0 :", scen_properties.R0_km)

    # Calculate the Gaussian distribution values for each shell
    gauss_values = np.array([np.exp(-(alt - mu_lam) ** 2 / sig_lam) for alt in scen_properties.R0_km]) # satellites/year

    # Normalize the Gaussian distribution to ensure total launch rate remains constant
    gauss_values_sum = np.sum(gauss_values)
    normalized_gauss_values = gauss_values / gauss_values_sum

    # Calculate the normalized launch rate for each shell
    for i, alt in enumerate(scen_properties.R0_km):
        h_inds = np.argmin(np.abs(h - alt))
        launch_rate_per_year = lam0 * normalized_gauss_values[i]
        launch_rate_per_time_step = launch_rate_per_year #/ time_steps_per_year
        launch_rate_matrix[h_inds, :] = launch_rate_per_time_step

    return launch_rate_matrix

def launch_sensitivity(h, scen_properties, mu_lam, sig_lam, lam0):
    """
    This function will return the lambda function for a required species.

    :param t: The time from the scenario start in years
    :type t: int
    :param h: The altitude above the ellipsoid in km of shell lower edge
    :type h: int
    :param species_properties: Species properties
    :type species_properties: Species
    :param scen_properties: Scenario properties
    :type scen_properties: ScenarioProperties
    :return: Lambdadot is the rate of change in the species in each sheel at the specified time due to launch
    :rtype: SciPy interp1d function
    """
    # # Find the index for the given altitude
    # h_inds = np.where(scen_properties.HMid == h)
    # print(species_properties.sym_name)

    # # Retrieve the appropriate lambda function for the altitude and evaluate it at time t
    # Lambdadot = species_properties.lambda_funs
    # return Lambdadot
    if len(h) != scen_properties.n_shells:
        raise ValueError("Constant launch rate must be specified per altitude shell.")

    n_shells = scen_properties.n_shells
    scen_times = scen_properties.scen_times
    simulation_duration = scen_properties.simulation_duration
    # Calculate the number of time steps per year
    time_steps_per_year = len(scen_times) / simulation_duration

    # Create a matrix filled with zeros
    launch_rate_matrix = np.zeros((n_shells, len(scen_times)))

    # Find the index of the shell closest to the specified altitude
    for alt in scen_properties.R0_km:

        h_inds = np.argmin(np.abs(h - alt))

        if h_inds < 0:
            # Calculate the number of satellites per time step
            launch_rate_per_year = 0  # satellites/year
            launch_rate_per_time_step = launch_rate_per_year / time_steps_per_year

            # Fill the row corresponding to the nearest altitude with the launch rate per time step
            launch_rate_matrix[h_inds, :] = launch_rate_per_time_step

        else:
            # Calculate the number of satellites per time step
            launch_rate_per_year = lam0*np.exp(-(alt-mu_lam)**2/sig_lam)  # satellites/year
            launch_rate_per_time_step = launch_rate_per_year / time_steps_per_year

            # Fill the row corresponding to the nearest altitude with the launch rate per time step
            launch_rate_matrix[h_inds, :] = launch_rate_per_time_step.round(4)

    return launch_rate_matrix


def julian_to_datetime(julian_date):
    # Julian Date for Unix epoch (1970-01-01)
    JULIAN_EPOCH = 2440587.5
    try:
        # Calculate the number of days from the Unix epoch
        days_from_epoch = julian_date - JULIAN_EPOCH
        # Create a datetime object for the Unix epoch and add the calculated days
        unix_epoch = datetime(1970, 1, 1)
        result_date = unix_epoch + timedelta(days=days_from_epoch)
        return result_date
    except OverflowError as e:
        # Handle dates that are out of range
        print(f"Date conversion error: {e}")
        return None

def BASE_init_model(scen_properties, file_path=None):
    """
    This function will create for the starting population

    :param scen_properties: Scneario properties
    :type scen_properties: ScenarioProperties
    :return: The initial population
    :rtype:  pandas.DataFrame
    """
    if file_path == None:
        x0 = pd.DataFrame(index=range(scen_properties.n_shells), columns=scen_properties.species_names)
        x0_params = [[650, 200 ** 2, 4000], [550, 250 ** 2, 800], [650, 200 ** 2, 400]]
        x0_params = pd.DataFrame(x0_params, index=scen_properties.species_names, columns=["mu", "sig", "O"])
        h = scen_properties.HMid
        # for species_name in scen_properties.species_names:
        #     print(species_name)
        #     print(x0_params["mu"][species_name])
        #     for alt in scen_properties.R0_km:
        #
        #         h_inds = np.argmin(np.abs(h - alt))
        #
        #         if h_inds < 0:
        #             # Calculate the number of satellites per time step
        #             init_population = 0  # initial population
        #
        #             # Fill the row corresponding to the nearest altitude with the launch rate per time step
        #             x0[species_name][h_inds] = init_population
        #
        #         else:
        #             # Calculate the number of satellites per time step
        #             init_population = x0_params["O"][species_name] * np.exp(
        #                 -(alt - x0_params["mu"][species_name]) ** 2 / x0_params["sig"][species_name])  # satellites/year
        #
        #             # Fill the row corresponding to the nearest altitude with the launch rate per time step
        #             x0[species_name][h_inds] = init_population


        for species_name in scen_properties.species_names:
            mult_factor = 1.6*(scipy.stats.norm.cdf(max(scen_properties.R0_km), loc=x0_params["mu"][species_name],
                                                scale=np.sqrt(x0_params["sig"][species_name])) - scipy.stats.norm.cdf(
                min(scen_properties.R0_km), loc=x0_params["mu"][species_name], scale=np.sqrt(x0_params["sig"][species_name])))
            init_gauss = np.array(
                [np.exp(-(alt - x0_params["mu"][species_name]) ** 2 / x0_params["sig"][species_name])
                 for alt in scen_properties.R0_km])
            init_gauss_sum = np.sum(init_gauss)
            init_gauss_norm = init_gauss / init_gauss_sum * mult_factor

            for i, alt in enumerate(scen_properties.R0_km):
                h_inds = np.argmin(np.abs(h - alt))
                init_population = x0_params["O"][species_name] * init_gauss_norm[h_inds]
                x0.at[h_inds, species_name] = init_population


        # for species_name in scen_properties.species_names:
        #     print(species_name)
        #     print(x0_params["mu"][species_name])

    else:
        x0 = pd.read_csv(file_path, names=scen_properties.species_names)
        #print(x0)

    return x0

def ADEPT_traffic_model(scen_properties, file_path):
    """
    From an initial population and future model csv, this function will create for the starting population, 
    then one for each time step in the future model.

    The output matrices will be in the form of a matrix, with the species as columns and the number of orbital shells as rows based on alt_bin.
    e.g. if you have 5 species and 10 shells, the matrix will be 10x5.

    :param scen_properties: Scneario properties
    :type scen_properties: ScenarioProperties
    :param file_path: Local File Path of CSV
    :type file_path: str
    :return: The initial population and future launch model
    :rtype:  pandas.DataFrame, pandas.DataFrame
    """
    # Load the traffic model data

    T = pd.read_csv(file_path)
    
    T['epoch_start_datime'] = T['epoch_start'].apply(lambda x: julian_to_datetime(x))

    if 'obj_class' not in T.columns:
        T = define_object_class(T)  # Make sure this function is defined and imported

    # Calculate Apogee, Perigee, and Altitude
    T['apogee'] = T['sma'] * (1 + T['ecc'])
    T['perigee'] = T['sma'] * (1 - T['ecc'])
    T['alt'] = (T['apogee'] + T['perigee']) / 2 - scen_properties.re

    # Map species type based on object class
    # species_dict = {"Non-station-keeping Satellite": "Sns",
    #                 "Rocket Body": "B",
    #                 "Station-keeping Satellite": "Su",
    #                 "Coordinated Satellite": "S",
    #                 "Debris": "N",
    #                 "Candidate Satellite": "C"}
    species_dict = {"Non-station-keeping Satellite": "S",
                    "Rocket Body": "N_223kg",
                    "Station-keeping Satellite": "S",
                    "Coordinated Satellite": "S",
                    "Debris": "N_0.64kg",
                    "Candidate Satellite": "C"}

    T['species_class'] = T['obj_class'].map(species_dict)

    #T = T[((T.species_class == "S") & (T.mass <= 600)) | ((T.species_class == "N_0.64kg") & (T.mass <= 1))| ((T.species_class == "N_223kg") & (T.mass <= 600))]

    # Initialize an empty DataFrame for new data
    T_new = pd.DataFrame()

    # Loop through object classes and assign species based on mass
    for obj_class in T['obj_class'].unique():
            species_class = species_dict.get(obj_class)
            
            if species_class in scen_properties.species_names:
                T_obj_class = T[T['obj_class'] == obj_class].copy()
                T_obj_class['species'] = species_class
                T_new = pd.concat([T_new, T_obj_class])

    # Assign objects to corresponding altitude bins
    T_new['alt_bin'] = T_new['alt'].apply(find_alt_bin, args=(scen_properties,))


    # Filter T_new to include only species present in scen_properties
    T_new = T_new[T_new['species_class'].isin(scen_properties.species_names)]

    # Initial population
    x0 = T_new[T_new['epoch_start_datime'] < scen_properties.start_date]

    # Create a pivot table, keep alt_bin
    df = x0.pivot_table(index='alt_bin', columns='species', aggfunc='size', fill_value=0)

    # Create a new data frame with column names like scenario_properties.species_sym_names and rows of length n_shells
    x0_summary = pd.DataFrame(index=range(scen_properties.n_shells), columns=scen_properties.species_names).fillna(0)
    x0_summary.index.name = 'alt_bin'

    # Merge the two dataframes
    for column in df.columns:
        if column in x0_summary.columns:
            x0_summary[column] = df[column]

    # fill NaN with 0
    x0_summary.fillna(0, inplace=True)

    # Future Launch Model
    flm_steps = pd.DataFrame()

    time_increment_per_step = scen_properties.simulation_duration / scen_properties.steps

    time_steps = [scen_properties.start_date + timedelta(days=365.25 * time_increment_per_step * i) 
                for i in range(scen_properties.steps + 1)]    

    for i, (start, end) in tqdm(enumerate(zip(time_steps[:-1], time_steps[1:])), total=len(time_steps)-1, desc="Processing Time Steps"):
        flm_step = T_new[(T_new['epoch_start_datime'] >= start) & (T_new['epoch_start_datime'] < end)]
        # print(f"Step: {start} - {end}, Objects: {flm_step.shape[0]}")
        flm_summary = flm_step.groupby(['alt_bin', 'species']).size().unstack(fill_value=0)

        # all objects aren't always in shells, so you need to these back in. 
        flm_summary = flm_summary.reindex(range(0, scen_properties.n_shells), fill_value=0)

        flm_summary.reset_index(inplace=True)
        flm_summary.rename(columns={'index': 'alt_bin'}, inplace=True)

        flm_summary['epoch_start_date'] = start # Add the start date to the table for reference
        flm_steps = pd.concat([flm_steps, flm_summary])

    #print("flm_steps: ", flm_steps)

    return x0_summary, flm_steps


def find_mass_bin(mass, scen_properties, species_cell):
    """
    Find the mass bin for a given mass.

    :param mass: Mass of the object in kg
    :type mass: float
    :param scen_properties: The scenario properties object
    :type scen_properties: ScenarioProperties
    :param species_cell: The species cell to find the mass bin for
    :type species_cell: Species
    :return: The mass bin for the given mass
    :rtype: int
    """
    for species in species_cell:
        if species.mass_lb <= mass < species.mass_ub:
            return species.sym_name
    return None

def find_alt_bin(altitude, scen_properties):
    """
    Given an altidude and the generic pySSEM properties, it will calculate the index from the R02 array

    :param altitude: Altitude of an object
    :type altitude: int
    :param scen_properties: The scenario properties object
    :type scen_properties: ScenarioProperties
    :return: Orbital Shell Array Index or None if out of range
    :rtype: int
    """
    shell_altitudes = scen_properties.R0_km

    # The case for an object where it is below the lowest altitude
    if altitude < shell_altitudes[0]:
        return
    
    # The case for an object where it is above the highest altitude
    if altitude >= shell_altitudes[-1]:
        return 

    for i in range(len(shell_altitudes)):  # -1 to prevent index out of range
        try:
            if shell_altitudes[i] <= altitude < shell_altitudes[i + 1]:
                return i  
        except IndexError: # This is the top most shell and will be the last one
            return len(shell_altitudes) 
    

def define_object_class(T):
    """
    Define the object class of each object in the traffic model.
    Adds them to a new column named "obj_type" or overwrites the existing column.

    :param T: list of launches
    :type T: pandas.DataFrame
    """

    T['obj_class'] = "Debris"

    # Classify Rocket Bodies
    T.loc[T['obj_type'] == 1, 'obj_class'] = "Rocket Body"

    # Classify Satellites
    T.loc[(T['obj_type'] == 2) & (T['stationkeeping'] != 0) & (T['stationkeeping'] < 5), 'obj_class'] = "Station-keeping Satellite"
    T.loc[(T['obj_type'] == 2) & (T['stationkeeping'] == 0), 'obj_class'] = "Non-station-keeping Satellite"
    T.loc[(T['obj_type'] == 2) & (T['stationkeeping'] == 5), 'obj_class'] = "Coordinated Satellite"
    T.loc[(T['obj_type'] == 2) & (T['stationkeeping'] == 6), 'obj_class'] = "Candidate Satellite"

    # Classify Debris
    T.loc[T['obj_type'] == 3 , 'obj_class'] = "Debris"

    # Count unclassified rows
    unclassed_rows = (T['obj_class'] == "Unknown").sum()
    if unclassed_rows > 0:
        print(f'\t{unclassed_rows} Unclassified rows remain.')

    return T