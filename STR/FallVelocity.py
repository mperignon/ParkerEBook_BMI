#! /usr/bin/env python
"""Basic Model Interface implementation for Parker's FallVelocity."""

import types

import numpy as np
import yaml
from basic_modeling_interface import Bmi

from STR_scripts.s_FallVelocity import FallVelocity_solver


class FallVelocity(Bmi):


    _name = 'STR_FallVelocity'
    
    # Standard Names of variables it could potentially get from another model
    _input_var_names = (
            'particle__diameter',
            'fluid__kinematic_viscosity',
            'gravitational_acceleration',
            'fluid__density',
            'particle__density',)
            
    # Standard Names of variables it could pass to another model
    _output_var_names = (
            'particle__settling_velocity',
            'particle__Reynolds_number',
            'particle__dimensionless_fall_velocity',)

    def __init__(self):
        """Create a model that is ready for initialization."""
        self._FV = None
        self._time = 0.
        self._values = {}
        self._var_units = {}



    def initialize(self, filename='../input_files/FallVelocity.yaml'):
        """Initialize the FallVelocity model.

        Parameters
        ----------
        filename : str, optional
            Path to name of input file.
        """
        
        with open(filename, 'r') as file_obj:
            params = yaml.load(file_obj)
        
        default_params = {'particle__diameter': 0.1,
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
                    
                    
            
        self._FV = FallVelocity_solver(params)


        self._values = {'particle__diameter': self._FV.grain_size,
                        'fluid__kinematic_viscosity': self._FV.kinematic_viscosity,
                        'gravitational_acceleration': self._FV.gravitational_acceleration,
                        'fluid__density': self._FV.density_of_fluid,
                        'particle__density': self._FV.density_of_particle,
                        'particle__settling_velocity': self._FV.settling_velocity,
                        'particle__Reynolds_number': self._FV.Reynolds_number,
                        'particle__dimensionless_fall_velocity': self._FV.dimensionless_fall_velocity}
        
        
        self._var_units = {
          'particle__diameter': 'mm',
          'fluid__kinematic_viscosity': 'm2 s-1',
          'gravitational_acceleration': 'm s-2',
          'fluid__density': 'Kg m-3',
          'particle__density': 'Kg m-3',
          'particle__settling_velocity': 'm s-1',
          'particle__Reynolds_number': '-',
          'particle__dimensionless_fall_velocity': '-'}
    


    # keep these methods so it can couple with other BMIed models
    def update(self):
        """Advance model by one time step.
           Since this script does not evolve in time,
           update() just runs the script once"""
           
        self._FV.run()

    def update_frac(self, time_frac):
        """Update model by a fraction of a time step.

        Parameters
        ----------
        time_frac : float
            Fraction fo a time step.
        
        Since this script does not evolve in time,
        update() just runs the script once"""
        
        self.update()


    def update_until(self, then):
        """Update model until a particular time.

        Parameters
        ----------
        then : float
            Time to run model until.
        
        Since this script does not evolve in time,
        update() just runs the script once"""
        
        self.update()

    def finalize(self):
        """Finalize model."""
        self._FV.finalize()
        self._FV = None
        
        
        



    def get_var_type(self, var_name):
        """Data type of variable.

        Parameters
        ----------
        var_name : str
            Name of variable as CSDMS Standard Name.

        Returns
        -------
        str
            Data type.
        """
        return str(self.get_value_ref(var_name).dtype)

    def get_var_units(self, var_name):
        """Get units of variable.

        Parameters
        ----------
        var_name : str
            Name of variable as CSDMS Standard Name.

        Returns
        -------
        str
            Variable units.
        """
        return self._var_units[var_name]

    def get_var_nbytes(self, var_name):
        """Get units of variable.

        Parameters
        ----------
        var_name : str
            Name of variable as CSDMS Standard Name.

        Returns
        -------
        int
            Size of data array in bytes.
        """
        return self.get_value_ref(var_name).nbytes



    def get_value_ref(self, var_name):
        """Reference to values.

        Parameters
        ----------
        var_name : str
            Name of variable as CSDMS Standard Name.

        Returns
        -------
        array_like
            Value array.
        """
        return self._values[var_name]

    def get_value(self, var_name):
        """Copy of values.

        Parameters
        ----------
        var_name : str
            Name of variable as CSDMS Standard Name.

        Returns
        -------
        array_like
            Copy of values.
        """
        return self.get_value_ref(var_name)


    def set_value(self, var_name, src):
        """Set model values.

        Parameters
        ----------
        var_name : str
            Name of variable as CSDMS Standard Name.
        src : int or float
            New value
        """
        self._values[var_name] = float(src)


    def get_component_name(self):
        """Name of the component."""
        return self._name

    def get_input_var_names(self):
        """Get names of input variables."""
        return self._input_var_names

    def get_output_var_names(self):
        """Get names of output variables."""
        return self._output_var_names
      
