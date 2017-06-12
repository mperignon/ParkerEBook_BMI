from math import sqrt, log10, pow
import json, os, yaml
import numpy as np
        

class FallVelocity_solver(object):
    """
    FallVelocity

    Python version of Gary Parker's 1D Sediment Transport Morphodynamics e-book,
    originally in Visual Basic and converted to C by Andrew Leman (2009)

    MC Perignon
    Oct 2015

    ----------------------------------------------

    Computes the settling velocity, vs, of a particle with the formulation of Dietrich (1982).

    Only valid for Reynold's numbers less than or equal to 2.5 x 10**6. If Re is greater than
    this upper limit, warns and exits.

    Dietrich, E. W., 1982, Settling velocity of natural particles,
    Water Resources Research, 18(6), 1626-1982.

    ----------------------------------------------
    Input:
    D: Sediment particle size (mm)
    nu: Kinematic viscosity of the liquid (m**2/s)
    g: Acceleration due to gravity (m/s**2)
    rho_w: Density of liquid (Kg/m**3)
    rho_s: Density of sediment(Kg/m**3)

    Output:
    vs: Particle settling velocity (m/s)
    Re: Particle Reynolds number
    Rf: Dimensionless fall velocity of the particle


    Equations:
    R = (rho_s - rho_w) / rho_w

    Re = sqrt(R*g*D)*D / nu

    Rf = vs / sqrt(R*g*D)
    """

    
    def __init__(self, params = None):

        if params is None:
        
            filename='input_files/FallVelocity.yaml'

            with open(filename, 'r') as file_obj:
                params = yaml.load(file_obj)

            default_params = {
                  'particle__diameter': 0.1,
                  'fluid__kinematic_viscosity': 1e-6,
                  'gravitational_acceleration': 9.81,
                  'fluid__density': 1000.,
                  'particle__density': 2650.,
                  'output_filename': 'output/FallVelocity_output.json'}


            # sets missing parameters to default values
            for key,value in default_params.items():
                params[key] = params.get(key, default_params[key])

            for key,value in params.items():
                if (value is None) or (value == 'None'):
                    params[key] = default_params[key]
        
        
        self._grain_size = np.array([float(params['particle_diameter'])]) # mm
        self._nu = np.array([float(params['fluid__kinematic_viscosity'])])
        self._g = np.array([float(params['gravitational_acceleration'])])
        self._rho_w = np.array([float(params['fluid__density'])])
        self._rho_s = np.array([float(params['particle__density'])])
        self._output_filename = str(params['output_filename'])
        
        # Outputs
        self._vs = np.empty((1,))
        self._Re = np.empty((1,))
        self._Rf = np.empty((1,))
    
    
    @property
    def grain_size(self):
        """Grain size in mm."""
        
        return self._grain_size


    @grain_size.setter
    def grain_size(self, new_D):
        """Set the grain size.

        Parameters
        ----------
        new_D : float
            New grain size in mm.
        """
        self._grain_size[:] = new_D
    
    @property
    def kinematic_viscosity(self):
        """Kinematic viscosity of the fluid in m2/s"""
        return self._nu
        
    @kinematic_viscosity.setter
    def kinematic_viscosity(self, new_nu):
        """Set the kinematic viscosity.

        Parameters
        ----------
        new_nu : float
            New kinematic viscosity in m2/s.
        """
        
        self._nu[:] = new_nu
        
        
    @property
    def gravitational_acceleration(self):
        """Gravitational acceleration in m/s2"""
        return self._g
        
    @gravitational_acceleration.setter
    def gravitational_acceleration(self, new_g):
        """Set the gravitational acceleration.

        Parameters
        ----------
        new_g : float
            New gravitational acceleration in m/s2.
        """
        
        self._g[:] = new_g
        
    @property
    def density_of_fluid(self):
        """Density of the fluid in Kg/m3"""
        return self._rho_w
        
    @density_of_fluid.setter
    def density_of_fluid(self, new_rhow):
        """Set the density of the fluid.

        Parameters
        ----------
        new_rhow : float
            New density of the fluid in Kg/m3.
        """
        
        self._rho_w[:] = new_rhow
        
    @property
    def density_of_particle(self):
        """Density of the particle in Kg/m3"""
        return self._rho_s


    @density_of_particle.setter
    def density_of_particle(self, new_rhos):
        """Set the density of the particle.

        Parameters
        ----------
        new_rhos : float
            New density of the particle in Kg/m3.
        """
        
        self._rho_s[:] = new_rhos
        
            
    
    @property
    def settling_velocity(self):
        """Settling velocity in m/s"""
        
        return self._vs
        
    
    @property
    def Reynolds_number(self):
        """Reynolds number of the particle"""
        
        return self._Re
    
    @property
    def dimensionless_fall_velocity(self):
        """Dimensionless fall velocity of the particle"""
        
        return self._Rf
        
        
            
    def run(self):
    
        self._R = (self._rho_s - self._rho_w) / self._rho_w
        
        self._Re[:] = self.calculateReynoldsnumber()     
        self._Rf[:] = self.calculate_dimensionless_fall_velocity()                                       
        self._vs[:] = self.calculate_settling_velocity()
        
        
        
    def calculateReynoldsnumber(self):

        _Re = sqrt(self._R * self._g * (self._grain_size/1000.)) * \
                    (self._grain_size/1000.) / self._nu

        failMessage  = "This equation is only valid for Reynolds numbers "\
                        "below 2.6e+06. The calculated Reynolds number for "\
                        "this particle is %.2g" % _Re
        assert _Re <= 2.6e6, failMessage       
        
        return _Re
  
  
        
    def calculate_dimensionless_fall_velocity(self):
        
        x = log10(self._Re**2)
        y = (-3.76715) + (1.92944*x) - (0.09815*x*x) - \
            (0.00575*x*x*x) + (0.00056*x*x*x*x)
        
        _Rf = pow(pow(10,y) / self._Re, 1./3)

        return _Rf

        
    def calculate_settling_velocity(self):
        
        _vs = self._Rf * sqrt(self._R * self._g * (self._grain_size/1000.))
        
        return _vs
        

    def _assure_path_exists(self, path):
        """Helper function to create output directories"""
        
        dir = os.path.dirname(path)
        if not os.path.exists(dir):
                os.makedirs(dir)        
            
            
            
    def finalize(self):

        
        self._assure_path_exists(self._output_filename)

        output_dict = {
            'Grain_size_mm' : list(self._grain_size),
            'Kinematic_viscosity' : list(self._nu),
            'Gravitational_acceleration' : list(self._g),
            'Density_of_fluid' : list(self._rho_w),
            'Density_of_sediment' : list(self._rho_s),
            'Reynolds_number' : list(self._Re),
            'Dimensionless_fall_velocity' : list(self._Rf),
            'Settling_velocity' : list(self._vs)
        }
        
        with open(self._output_filename, 'w') as f:
            json.dump(output_dict, f, indent=4)

                
            
            
            
