#/////////////////////////////////////////////////////////////
#// Main author: Chengshun Shang (CIMNE)
#// Email: cshang@cimne.upc.edu; chengshun.shang1996@gmail.com
#// Date: July 2024
#/////////////////////////////////////////////////////////////

from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

import KratosMultiphysics.DEMApplication.sphere_strategy as SolverStrategy
BaseExplicitStrategy = SolverStrategy.ExplicitStrategy

import math

class ExplicitStrategy(BaseExplicitStrategy):

    def __init__(self, all_model_parts, creator_destructor, dem_fem_search, DEM_parameters, procedures):

        BaseExplicitStrategy.__init__(self, all_model_parts, creator_destructor, dem_fem_search, DEM_parameters, procedures)

        #if "PostSkinSphere" in DEM_parameters.keys():
        #    self.print_skin_sphere = DEM_parameters["PostSkinSphere"].GetBool()
        self.polyhedron_model_part = all_model_parts.Get("PolyhedronPart")


    def SetVariablesAndOptions(self):

        BaseExplicitStrategy.SetVariablesAndOptions(self)

        self.settings = ContactExplicitSolverSettings()
        self.settings.r_model_part = self.spheres_model_part
        self.settings.contact_model_part = self.contact_model_part
        self.settings.fem_model_part = self.fem_model_part
        self.settings.inlet_model_part = self.inlet_model_part
        self.settings.cluster_model_part = self.cluster_model_part
        self.settings.polyhedron_model_part = self.polyhedron_model_part

    def CreateCPlusPlusStrategy(self):

        self.SetVariablesAndOptions()

        # ADDITIONAL VARIABLES AND OPTIONS
        #self.spheres_model_part.ProcessInfo.SetValue(CONTINUUM_SEARCH_RADIUS_AMPLIFICATION_FACTOR, self.continuum_search_radius_amplification_factor)

        '''
        for properties in self.spheres_model_part.Properties:
            for subproperties in properties.GetSubProperties():
                if subproperties.Has(DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME):
                    continuum_constitutive_law_name = subproperties[DEM_CONTINUUM_CONSTITUTIVE_LAW_NAME]
                    continuum_constitutive_law_instance = globals().get(continuum_constitutive_law_name)()
                    if continuum_constitutive_law_instance.CheckRequirementsOfStressTensor():
                        self.spheres_model_part.ProcessInfo.SetValue(COMPUTE_STRESS_TENSOR_OPTION, 1)
                        break'''

        self.cplusplus_strategy = ContactExplicitSolverStrategy(self.settings, self.max_delta_time, self.n_step_search, self.safety_factor,
                                                  self.delta_option, self.creator_destructor, self.dem_fem_search, self.search_strategy, self.solver_settings)

    def BeforeInitialize(self):
        self.CreateCPlusPlusStrategy()
        self.RebuildListOfPolyhedronParticles()
        self.SetNormalRadiiOnAllParticles()
        self.SetSearchRadiiOnAllParticles()

    def Initialize(self):
        self.cplusplus_strategy.Initialize()  # Calls the cplusplus_strategy Initialize function (initializes all elements and performs other necessary tasks before starting the time loop) (C++)

    def AddAdditionalVariables(self, polyhedron_model_part, DEM_parameters):
        #spheres_model_part.AddNodalSolutionStepVariable(COHESIVE_GROUP) 
        pass

    def RebuildListOfPolyhedronParticles(self):
        self.cplusplus_strategy.RebuildListOfPolyhedronParticles()

    def SetNormalRadiiOnAllParticles(self):
        self.cplusplus_strategy.SetNormalRadiiOnAllParticles(self.polyhedron_model_part)
