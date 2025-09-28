
from ..launch.launch import ADEPT_traffic_model, launch_func_constant, BASE_init_model
from ..pmd.pmd import pmd_func_sat, pmd_func_derelict
from ..drag.drag import densityexp
from .scen_properties import ScenarioProperties
from .species_pair_class import SpeciesPairClass

import numpy as np
from sympy import zeros, symbols, sqrt, exp, pi


def define_launch_events(species_list, model):
    events = []

    launch_rate_matrix = species_list[0].launch_func(
        model.scenario_properties.scen_times,
        model.scenario_properties.HMid,
        species_list[0],
        model.scenario_properties
    )  # shape: (n_shells, n_time_steps), objects/time_step

    # calculate time_steps_per_year
    time_steps_per_year = len(model.scenario_properties.scen_times) #/ model.scenario_properties.simulation_duration
    print(time_steps_per_year)

    for shell_idx in range(model.scenario_properties.n_shells):
        
        print("launch_mtx: ", launch_rate_matrix[shell_idx, 0])

        # calculate yearly rates
        launch_rate_per_step = launch_rate_matrix[shell_idx, 0]  # first step (objects/time_step)
        launch_rate_per_year = launch_rate_per_step # convert to objects/year

        # total_launches_per_shell = np.sum(launch_rate_matrix, axis=1)
        # launch_rate_per_year = total_launches_per_shell / len(model.scenario_properties.scen_times) #/ model.scenario_properties.simulation_duration
        # print("launch_rate_post", launch_rate_per_year)

        events.append({
            'name': f'satellite_launch_shell_{shell_idx+1}',
            'rate': launch_rate_per_year,
            'jump': {f'shell_{shell_idx+1}_S': +1}
        })

    print("launch")
    for event in events:
        print(event)

    return events
    
# def define_launch_events(species_list, model):
#     events = []

#     launch_rate_matrix = species_list[0].launch_func(
#         model.scenario_properties.scen_times,
#         model.scenario_properties.HMid,
#         species_list[0],
#         model.scenario_properties
#     )

#     total_launches_per_shell = np.sum(launch_rate_matrix, axis=1)
#     simulation_years = model.scenario_properties.simulation_duration
#     launch_rate_per_year = total_launches_per_shell / simulation_years

#     for shell_idx in range(model.scenario_properties.n_shells):
#         events.append({
#             'name': f'satellite_launch_shell_{shell_idx+1}',
#             'rate': launch_rate_per_year[shell_idx],
#             'jump': {f'shell_{shell_idx+1}_S': +1}
#         })

#     return events

def define_pmd_events(species_properties, model, state):
    events = []

    S_deltat = species_properties[0].deltat  # Satellite lifetime
    shells = model.scenario_properties.n_shells
    
    for k in range(shells):
        S_pop = state[f'shell_{k+1}_S']  
        
        # PMD success event (Satellite decrease)
        pmd_S_rate = S_pop / S_deltat
        events.append({
            'name': f'pmd_sat_success_shell_{k+1}',
            'rate': pmd_S_rate,
            'jump': {f'shell_{k+1}_S': -1}
        })

        # PMD fail event (Satellite decrease & Derelict increase)
        for species in species_properties[2].pmd_linked_species:
            Pm = species.Pm  # PMD success prob 
            pmd_fail_rate = (1 - Pm) * (S_pop / S_deltat)
            events.append({
                'name': f'pmd_fail_increase_derelict_shell_{k+1}',
                'rate': pmd_fail_rate,
                'jump': {f'shell_{k+1}_S': -1, f'shell_{k+1}_N_223kg': +1}
            })
            
    # print("pmd")
    # for event in events:
    #     print(event)
        
    return events

def update_pmd_events(all_species, model, state):
    species_properties = all_species 
    return define_pmd_events(species_properties, model, state)

def define_collision_events_from_pairs(model, state):
    events = []
    
    n_f = symbols(f'n_f:{model.scenario_properties.n_shells}')
    
    for pair in model.scenario_properties.collision_pairs:
        phi_values = np.array(pair.phi).flatten()
        gamma_values = np.array(pair.gammas)

        species1_name = pair.species1.sym_name
        species2_name = pair.species2.sym_name

        for shell_idx in range(model.scenario_properties.n_shells):
            species1_pop = state[f'shell_{shell_idx+1}_{species1_name}']
            species2_pop = state[f'shell_{shell_idx+1}_{species2_name}']

            base_collision_rate = phi_values[shell_idx] * species1_pop * species2_pop

            gamma_shell = gamma_values[shell_idx]

            nf_val = np.array(pair.nf).flatten()

            # init event
            lethal_jump = {}
            disabling_jump = {}
            lethal_rate = 0
            disabling_rate = 0
            jump = {}
            rate = 0

            # debris numeric
            debris_jump = int(gamma_shell[2].subs(n_f[shell_idx], nf_val[shell_idx]))


            # Species 1 처리 (Satellite 또는 Derelict)
            if species1_name == 'S':
                lethal_jump[f'shell_{shell_idx+1}_{species1_name}'] = -1  # Satellite decrease
                disabling_jump[f'shell_{shell_idx+1}_{species1_name}'] = -1  # Satellite decrease

                if species2_name == 'S':
                    lethal_jump[f'shell_{shell_idx+1}_{species2_name}'] = -1  # Satellite decrease
                    lethal_rate = base_collision_rate * pair.species1.alpha_active

                    events.append({
                        'name': f'collision_{species1_name}_{species2_name}_lethal_shell_{shell_idx+1}',
                        'rate': lethal_rate,
                        'jump': {f'shell_{shell_idx+1}_{species1_name}': -2, f'shell_{shell_idx+1}_{species2_name}': -2, f'shell_{shell_idx+1}_N_0.64kg': debris_jump/pair.species1.alpha_active}
                    })

                elif species2_name == 'N_223kg':
                    lethal_jump[f'shell_{shell_idx+1}_{species2_name}'] = -1  # Derelict decrease
                    disabling_jump[f'shell_{shell_idx+1}_{species2_name}'] = 1  # Derelict increase
                    lethal_rate = base_collision_rate * pair.species1.alpha
                    disabling_rate = base_collision_rate * pair.species1.delta
                    
                    events.append({
                        'name': f'collision_{species1_name}_{species2_name}_lethal_shell_{shell_idx+1}',
                        'rate': lethal_rate,
                        'jump': {f'shell_{shell_idx+1}_{species1_name}': -1, f'shell_{shell_idx+1}_{species2_name}': -1, f'shell_{shell_idx+1}_N_0.64kg': debris_jump/pair.species1.alpha}
                    })

                    events.append({
                        'name': f'collision_{species1_name}_{species2_name}_disabling_shell_{shell_idx+1}',
                        'rate': disabling_rate,
                        'jump': {f'shell_{shell_idx+1}_{species1_name}': -1, f'shell_{shell_idx+1}_{species2_name}': 1}
                    })
                
                elif species2_name == 'N_0.64kg':
                    lethal_jump[f'shell_{shell_idx+1}_{species2_name}'] = gamma_shell[2].subs(n_f[shell_idx], nf_val[shell_idx])  # Debris increase
                    disabling_jump[f'shell_{shell_idx+1}_{species2_name}'] = 0  # Derelict increase
                    lethal_rate = base_collision_rate * pair.species1.alpha
                    disabling_rate = base_collision_rate * pair.species1.delta
                    
                    events.append({
                        'name': f'collision_{species1_name}_{species2_name}_lethal_shell_{shell_idx+1}',
                        'rate': lethal_rate,
                        'jump': {f'shell_{shell_idx+1}_{species1_name}': -1, f'shell_{shell_idx+1}_N_0.64kg': debris_jump/pair.species1.alpha}
                    })

                    events.append({
                        'name': f'collision_{species1_name}_{species2_name}_disabling_shell_{shell_idx+1}',
                        'rate': disabling_rate,
                        'jump': {f'shell_{shell_idx+1}_{species1_name}': -1, f'shell_{shell_idx+1}_N_223kg': 1}
                    })
                    
                    
            elif species1_name == 'N_223kg':
                jump_amount = int(gamma_shell[0])
                jump[f'shell_{shell_idx+1}_{species1_name}'] = jump_amount
                rate = base_collision_rate  # lethal derelict collision

                if species2_name == 'N_223kg':
                    jump[f'shell_{shell_idx+1}_{species2_name}'] = -1  # Derelict decrease
                    rate = base_collision_rate
                    
                    events.append({
                        'name': f'collision_{species1_name}_{species2_name}_shell_{shell_idx+1}',
                        'rate': rate,
                        'jump': {f'shell_{shell_idx+1}_{species1_name}': -2, f'shell_{shell_idx+1}_{species2_name}': -2, f'shell_{shell_idx+1}_N_0.64kg': debris_jump}
                    })
                    
                elif species2_name == 'N_0.64kg':
                    jump[f'shell_{shell_idx+1}_{species2_name}'] = gamma_shell[2].subs(n_f[shell_idx], nf_val[shell_idx])  # Debris increase\
                    rate = base_collision_rate
                    
                    events.append({
                        'name': f'collision_{species1_name}_{species2_name}_shell_{shell_idx+1}',
                        'rate': rate,
                        'jump': {f'shell_{shell_idx+1}_{species1_name}': -1, f'shell_{shell_idx+1}_N_0.64kg': debris_jump}
                    })

            elif species1_name == 'N_0.64kg':
                jump_amount = int(gamma_shell[0])
                jump[f'shell_{shell_idx+1}_{species1_name}'] = jump_amount
                rate = base_collision_rate  # lethal derelict collision

                if species2_name == 'N_223kg':
                    jump[f'shell_{shell_idx+1}_{species2_name}'] = -1  # Derelict decrease
                    rate = base_collision_rate
                    
                    events.append({
                        'name': f'collision_{species1_name}_{species2_name}_shell_{shell_idx+1}',
                        'rate': rate,
                        'jump': {f'shell_{shell_idx+1}_{species2_name}': -1, f'shell_{shell_idx+1}_N_0.64kg': debris_jump}
                    })
                    
                elif species2_name == 'N_0.64kg':
                    jump[f'shell_{shell_idx+1}_{species2_name}'] = gamma_shell[2].subs(n_f[shell_idx], nf_val[shell_idx])  # Debris increase\
                    rate = base_collision_rate
                    
                    events.append({
                        'name': f'collision_{species1_name}_{species2_name}_shell_{shell_idx+1}',
                        'rate': rate,
                        'jump': {f'shell_{shell_idx+1}_N_0.64kg': gamma_shell[2].subs(n_f[shell_idx], nf_val[shell_idx])}
                    })

    print("colls")
    for event in events:
        print(event)
    
    return events





def update_collision_event_rates(collision_pairs, state, n_shells):
    collision_events = []

    n_f = symbols(f'n_f:{n_shells}')
    
    for pair in collision_pairs:
        phi_values = np.array(pair.phi).flatten()  
        gamma_values = np.array(pair.gammas)

        species1_name = pair.species1.sym_name
        species2_name = pair.species2.sym_name

        for shell_idx in range(n_shells):
            species1_pop = state[f'shell_{shell_idx+1}_{species1_name}']
            species2_pop = state[f'shell_{shell_idx+1}_{species2_name}']

            base_collision_rate = phi_values[shell_idx] * species1_pop * species2_pop
            gamma_shell = gamma_values[shell_idx]

            nf_val = np.array(pair.nf).flatten()

            # Initialize jumps and rates
            lethal_jump = {}
            disabling_jump = {}
            lethal_rate = 0
            disabling_rate = 0
            jump = {}
            rate = 0

            # debris numeric
            debris_jump = int(gamma_shell[2].subs(n_f[shell_idx], nf_val[shell_idx]))


            # Species 1 processing
            if species1_name == 'S':
                lethal_jump[f'shell_{shell_idx+1}_{species1_name}'] = -1
                disabling_jump[f'shell_{shell_idx+1}_{species1_name}'] = -1

                if species2_name == 'S':
                    lethal_jump[f'shell_{shell_idx+1}_{species2_name}'] = -1
                    lethal_rate = base_collision_rate * pair.species1.alpha_active

                    collision_events.append({
                        'name': f'collision_{species1_name}_{species2_name}_lethal_shell_{shell_idx+1}',
                        'rate': lethal_rate,
                        'jump': {
                            f'shell_{shell_idx+1}_{species1_name}': -2,
                            f'shell_{shell_idx+1}_{species2_name}': -2,
                            f'shell_{shell_idx+1}_N_0.64kg': debris_jump/pair.species1.alpha_active
                        }
                    })

                elif species2_name == 'N_223kg':
                    lethal_rate = base_collision_rate * pair.species1.alpha
                    disabling_rate = base_collision_rate * pair.species1.delta

                    # Lethal event
                    collision_events.append({
                        'name': f'collision_{species1_name}_{species2_name}_lethal_shell_{shell_idx+1}',
                        'rate': lethal_rate,
                        'jump': {
                            f'shell_{shell_idx+1}_{species1_name}': -1,
                            f'shell_{shell_idx+1}_{species2_name}': -1,
                            f'shell_{shell_idx+1}_N_0.64kg': debris_jump/pair.species1.alpha
                        }
                    })

                    # Disabling event
                    collision_events.append({
                        'name': f'collision_{species1_name}_{species2_name}_disabling_shell_{shell_idx+1}',
                        'rate': disabling_rate,
                        'jump': {
                            f'shell_{shell_idx+1}_{species1_name}': -1,
                            f'shell_{shell_idx+1}_{species2_name}': 1
                        }
                    })

                elif species2_name == 'N_0.64kg':
                    lethal_rate = base_collision_rate * pair.species1.alpha
                    disabling_rate = base_collision_rate * pair.species1.delta

                    # Lethal event
                    collision_events.append({
                        'name': f'collision_{species1_name}_{species2_name}_lethal_shell_{shell_idx+1}',
                        'rate': lethal_rate,
                        'jump': {
                            f'shell_{shell_idx+1}_{species1_name}': -1,
                            f'shell_{shell_idx+1}_N_0.64kg': debris_jump/pair.species1.alpha
                        }
                    })

                    # Disabling event
                    collision_events.append({
                        'name': f'collision_{species1_name}_{species2_name}_disabling_shell_{shell_idx+1}',
                        'rate': disabling_rate,
                        'jump': {
                            f'shell_{shell_idx+1}_{species1_name}': -1,
                            f'shell_{shell_idx+1}_N_223kg': 1
                        }
                    })

            elif species1_name == 'N_223kg':
                if species2_name == 'N_223kg':
                    collision_events.append({
                        'name': f'collision_{species1_name}_{species2_name}_shell_{shell_idx+1}',
                        'rate': base_collision_rate,
                        'jump': {
                            f'shell_{shell_idx+1}_{species1_name}': -2,
                            f'shell_{shell_idx+1}_{species2_name}': -2,
                            f'shell_{shell_idx+1}_N_0.64kg': debris_jump
                        }
                    })
                elif species2_name == 'N_0.64kg':
                    collision_events.append({
                        'name': f'collision_{species1_name}_{species2_name}_shell_{shell_idx+1}',
                        'rate': base_collision_rate,
                        'jump': {
                            f'shell_{shell_idx+1}_{species1_name}': -1,
                            f'shell_{shell_idx+1}_N_0.64kg': debris_jump
                        }
                    })

            elif species1_name == 'N_0.64kg':
                if species2_name == 'N_223kg':
                    collision_events.append({
                        'name': f'collision_{species1_name}_{species2_name}_shell_{shell_idx+1}',
                        'rate': base_collision_rate,
                        'jump': {
                            f'shell_{shell_idx+1}_{species2_name}': -1,
                            f'shell_{shell_idx+1}_N_0.64kg': debris_jump
                        }
                    })
                elif species2_name == 'N_0.64kg':
                    collision_events.append({
                        'name': f'collision_{species1_name}_{species2_name}_shell_{shell_idx+1}',
                        'rate': base_collision_rate,
                        'jump': {
                            f'shell_{shell_idx+1}_N_0.64kg': debris_jump
                        }
                    })

    return collision_events



def calculate_drag_event_rates(species, scen_properties, state):
    seconds_per_year = 365.25 * 24 * 3600
    mu = scen_properties.mu
    re = scen_properties.re
    R0_km = scen_properties.R0_km
    deltaH = scen_properties.deltaH

    rho = densityexp(R0_km)

    upper_drag_rates = np.zeros(scen_properties.n_shells)
    current_drag_rates = np.zeros(scen_properties.n_shells)

    if species.drag_effected:
        for k in range(scen_properties.n_shells):
            species_population_current = state[f'shell_{k+1}_{species.sym_name}']
            
            rvel_current = float(species.beta * sqrt(scen_properties.mu * scen_properties.R0[k+1]) * seconds_per_year)
            current_drag_rates[k] = rvel_current * rho[k] * species_population_current / deltaH[k]

            # upper shell --> current shell 
            if k < scen_properties.n_shells - 1:
                species_population_upper = state[f'shell_{k+2}_{species.sym_name}']
                
                rvel_upper = float(species.beta * sqrt(scen_properties.mu * scen_properties.R0[k]) * seconds_per_year)
                upper_drag_rates[k] = rvel_upper * rho[k+1] * species_population_upper / deltaH[min(k+1, scen_properties.n_shells-1)]

    return upper_drag_rates, current_drag_rates


def define_drag_events(species_list, model, state):
    events = []

    for species in species_list:
        upper_drag_rates, current_drag_rates = calculate_drag_event_rates(species, model.scenario_properties, state)

        for k in range(model.scenario_properties.n_shells):

            # current shell --> lower shell
            if k > 0:
                events.append({
                    'name': f'drag_down_shell_{k+1}_to_{k}',
                    'rate': current_drag_rates[k],
                    'jump': {
                        f'shell_{k+1}_{species.sym_name}': -1,
                        f'shell_{k}_{species.sym_name}': +1
                    }
                })
            else:
                # objects in the lowest shell removed
                events.append({
                    'name': 'drag_out_shell_1_to_atmosphere',
                    'rate': current_drag_rates[k],
                    'jump': {
                        f'shell_{k+1}_{species.sym_name}': -1
                    }
                })

            # upper shell --> current shell 
            if k < model.scenario_properties.n_shells - 1:
                events.append({
                    'name': f'drag_from_upper_shell_{k+2}_to_{k+1}',
                    'rate': upper_drag_rates[k],
                    'jump': {
                        f'shell_{k+2}_{species.sym_name}': -1,
                        f'shell_{k+1}_{species.sym_name}': +1
                    }
                })

    print("drag")
    for event in events:
        print(event)
        
    return events

def update_drag_events(species_list, model, state):
    updated_events = []

    for species in species_list:
        upper_drag_rates, current_drag_rates = calculate_drag_event_rates(species, model.scenario_properties, state)

        for k in range(model.scenario_properties.n_shells):

            if k > 0:
                updated_events.append({
                    'name': f'drag_down_shell_{k+1}_to_{k}',
                    'rate': current_drag_rates[k],
                    'jump': {
                        f'shell_{k+1}_{species.sym_name}': -1,
                        f'shell_{k}_{species.sym_name}': +1
                    }
                })
            else:
                updated_events.append({
                    'name': 'drag_out_shell_1_to_atmosphere',
                    'rate': current_drag_rates[k],
                    'jump': {
                        f'shell_{k+1}_{species.sym_name}': -1
                    }
                })

            if k < model.scenario_properties.n_shells - 1:
                updated_events.append({
                    'name': f'drag_from_upper_shell_{k+2}_to_{k+1}',
                    'rate': upper_drag_rates[k],
                    'jump': {
                        f'shell_{k+2}_{species.sym_name}': -1,
                        f'shell_{k+1}_{species.sym_name}': +1
                    }
                })

    return updated_events
