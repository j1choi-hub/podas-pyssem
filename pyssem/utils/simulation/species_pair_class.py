from sympy import symbols, Matrix, pi, S, Expr, zeros, sqrt
import matplotlib.pyplot as plt
import numpy as np

class SpeciesPairClass:
    def __init__(self, species1, species2, gammas, gammas_sqrt, source_sinks, scen_properties, model_type, fragsMadeDV=None):
        """
        This makes the species pair class associated with a collision between species1
        and species2. It will then create equations for the collision probability modifiers
        in gamma and the species in source_sinks.

        If the symbolic argument "n_f" is passed, it will be replaced with the n_f value
        for a collision involved species1 and species2 at each dv in scen_properties.v_imp2

        Other species will gain from a collision - e.g. debris. 

        Args:
            species1 (Species): The first species in the collision
            species2 (Species): The second species in the collision
            gammas (np.ndarray): The collision probability modifiers for each species in source_sinks.
            A scalar or a N x M matrix, where N is the number of altitude bins and M is the number of species
            with population addition/subtractions where this collision types occur. 
            source_sinks (list): A list of species that are either sources or sinks in the collision
            scen_properties (ScenarioProperties): The scenario properties object
        """

        if gammas.shape[1] != len(source_sinks):
            raise ValueError("Gammas and source_sinks must be the same length")

        # As species is a dictionary, it needs to be flatted first        
        all_species = [species for category in scen_properties.species.values() for species in category]

        time_steps_per_year = scen_properties.steps / scen_properties.simulation_duration

        self.name = f"species_pair({species1.sym_name}, {species2.sym_name})"
        print("pairing", {species1.sym_name}, {species2.sym_name})
        self.species1 = species1
        self.species2 = species2

        meter_to_km = 1 / 1000

        # Square of impact parameter
        if model_type == "JASON":
            self.sigma = 0.01
        else:
            self.sigma = (species1.radius + species2.radius) ** 2

        # Scaling based on v_imp, shell volume, and object radii
        self.phi = 3.1415 * scen_properties.v_imp2 / (scen_properties.V) * self.sigma * S(86400) * S(365.25) #/time_steps_per_year

        # Check if collision is catastrophic
        self.catastrophic = self.is_catastrophic(species1, species2, scen_properties.v_imp2)

        # Fragment generation equations
        M1 = species1.mass
        M2 = species2.mass
        LC = scen_properties.LC

        nf = zeros(len(scen_properties.v_imp2), 1)

        for i, dv in enumerate(scen_properties.v_imp2):
            if self.catastrophic[i]:
                # number of fragments generated during a catastrophic collision (NASA standard break-up model). M is the sum of the mass of the objects colliding in kg
                n_f_catastrophic = 0.1 * LC**(-S(1.71)) * (M1 + M2)**(S(0.75))
                nf[i] = n_f_catastrophic
            else:
                # number of fragments generated during a non-catastrophic collision (improved NASA standard break-up model: takes into account the kinetic energy). M is the mass of the less massive object colliding in kg
                n_f_damaging = 0.1 * LC**(-S(1.71)) * (min(M1, M2) * dv**2)**(S(0.75))
                nf[i] = n_f_damaging

        self.nf = nf.transpose()
            
        self.gammas = gammas
        self.gammas_sqrt = gammas_sqrt
        self.source_sinks = source_sinks
        self.eqs = Matrix(scen_properties.n_shells, len(all_species), lambda i, j: 0)
        self.eqs_diffusion = Matrix(scen_properties.n_shells, len(all_species), lambda i, j: 0)

        if isinstance(self.phi, (int, float, Expr)):
            phi_matrix = Matrix([self.phi] * len(gamma))
        else:
            phi_matrix = Matrix(self.phi)

        if scen_properties.fragment_spreading:
            product_sym = species1.sym.multiply_elementwise(species2.sym).T


        # Go through each gamma (which modifies collision for things like collision avoidance, or fragmentation into 
        # derelicsts, etc.) We increment the eqs matrix with the gamma * phi * species1 * species2.
        for i in range(gammas.shape[1]):
            gamma = gammas[:, i]
            gamma_sqrt = gammas_sqrt[:, i]
            #print(source_sinks[i].sym_name, gamma)
            eq_index = None
            for idx, spec in enumerate(all_species):
                if spec.sym_name == source_sinks[i].sym_name:
                    eq_index = idx
                    break

            if eq_index is None:
                raise ValueError(f"Equation index not found for {source_sinks[i].sym_name}")

            n_f = symbols(f'n_f:{scen_properties.n_shells}')

            if scen_properties.fragment_spreading:
                if i < 2:  # As first two columns are the reduction of the species in the collision (i.e -1)
                    eq = gamma.multiply_elementwise(phi_matrix).multiply_elementwise(species1.sym).multiply_elementwise(species2.sym)
                    eq_diffusion = gamma_sqrt.multiply_elementwise(phi_matrix.applyfunc(sqrt)).multiply_elementwise(species1.sym.applyfunc(sqrt)).multiply_elementwise(species2.sym.applyfunc(sqrt))
                else:  # Debris generated from collision
                    try:                         
                        fragsMadeDVcurrentDeb = fragsMadeDV[:, i-2] # First two rows are the reduction of the species in the collision (i.e -1)

                        # Create the 2D fragment matrix with circular shifts
                        fragsMade2D_list = [np.roll(fragsMadeDVcurrentDeb, shift) for shift in range(scen_properties.n_shells + 1)]
                        fragsMade2D = np.column_stack(fragsMade2D_list)

                        # Adjust the slicing to match MATLAB's slicing
                        fragsMade2D = fragsMade2D[scen_properties.n_shells:, :scen_properties.n_shells]  # from N_shell:end for rows, 1:N_shell for columns
                        fragsMade2D_sym = Matrix(fragsMade2D)

                        # Use the species product matrix and repeat it for each shell, this will allow all symbolic variables to be affected across all shells
                        rep_mat_sym = Matrix.vstack(*[product_sym for _ in range(scen_properties.n_shells)])

                        # Perform element-wise multiplication
                        sum_ = fragsMade2D_sym.multiply_elementwise(rep_mat_sym)
                                          
                        # Sum the columns of the multiplied_matrix
                        sum_matrix = Matrix([sum(sum_[row, :]) for row in range(sum_.shape[0])])

                        # Multiply gammas, phi, and the sum_matrix element-wise
                        eq = -gammas[:, 0].multiply_elementwise(phi_matrix).multiply_elementwise(sum_matrix)
                        eq_diffusion = -gammas_sqrt[:, 0].multiply_elementwise(phi_matrix.applyfunc(sqrt)).multiply_elementwise(sum_matrix.applyfunc(sqrt))

                        # Plotting (similar to MATLAB's imagesc)
                        # plt.figure(100)
                        # plt.clf()
                        # plt.imshow(fragsMade2D, aspect='auto', interpolation='none')
                        # plt.colorbar()
                        # plt.title(f"{self.name} for {source_sinks[i].sym_name}", fontsize=10)
                        # plt.gca().invert_yaxis()
                        # plt.show()
                        # print(f"eq: {eq}")
                    except Exception as e:
                        if fragsMadeDV == 0:
                            eq = gamma.multiply_elementwise(phi_matrix).multiply_elementwise(species1.sym).multiply_elementwise(species2.sym)
                            eq_diffusion = gamma_sqrt.multiply_elementwise(phi_matrix.applyfunc(sqrt)).multiply_elementwise(species1.sym.applyfunc(sqrt)).multiply_elementwise(species2.sym.applyfunc(sqrt))
                            continue
                        print(f"Error in creating debris matrix: {e}")
            else:
                eq = gamma.multiply_elementwise(phi_matrix).multiply_elementwise(species1.sym).multiply_elementwise(species2.sym)
                eq_diffusion = gamma_sqrt.multiply_elementwise(phi_matrix.applyfunc(sqrt)).multiply_elementwise(species1.sym.applyfunc(sqrt)).multiply_elementwise(species2.sym.applyfunc(sqrt))
            for j, val in enumerate(self.nf):
                eq = eq.subs(n_f[j], val)
                eq_diffusion = eq_diffusion.subs(n_f[j], val)

            self.eqs[:, eq_index] = self.eqs[:, eq_index] + eq
            self.eqs_diffusion[:, eq_index] = self.eqs_diffusion[:, eq_index] + eq_diffusion
            # print("eqs for :", species1.sym_name, species2.sym_name, self.eqs)
            # print("eqd for :", species1.sym_name, species2.sym_name, self.eqs_diffusion)



    def is_catastrophic(self, species1, species2, vels):
        """
        Determines if a collision is catastropic or non-catastrophic by calculating the 
        relative kinetic energy. If the energy is greater than 40 J/g, the collision is
        catastrophic from Johnson et al. 2001 (Collision Section)

        Args:
            mass1 (float): mass of species 1, kg
            mass2 (float): mass of species 2, kg
            vels (np.ndarray): array of the relative velocities (km/s) for each shell
        
        Returns:
            shell-wise list of bools (true if catastrophic, false if not catastrophic)
        """

        mass1 = species1.mass
        mass2 = species2.mass

        if mass1 <= mass2:
            smaller_mass = mass1
            larger_mass = mass2
        else:
            smaller_mass = mass2
            larger_mass = mass1

        
        # smaller_mass_g = smaller_mass * (1000) # kg to g

        # energy = [0.5 * smaller_mass * (v)**2 for v in vels] # Need to also convert km/s to m/s
        # is_catastrophic = [True if e/smaller_mass_g > 40 else False for e in energy]
        #
        # collision_mass = [smaller_mass for v in vels] # list of collision mass
        # is_catastrophic = [True if c_mass > 100 else False for c_mass in collision_mass]

        collision_mass = [smaller_mass for v in vels] # list of collision mass
        is_catastrophic = [True if c_mass > 100 else False for c_mass in collision_mass]

        return is_catastrophic
