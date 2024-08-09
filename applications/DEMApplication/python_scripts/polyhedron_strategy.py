#/////////////////////////////////////////////////////////////
#// Main author: Chengshun Shang (CIMNE)
#// Email: cshang@cimne.upc.edu; chengshun.shang1996@gmail.com
#// Date: July 2024
#/////////////////////////////////////////////////////////////

from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

import KratosMultiphysics.DEMApplication.sphere_strategy as SolverStrategy
BaseExplicitStrategy = SolverStrategy.ExplicitStrategy

import KratosMultiphysics.DEMApplication.polyhedron_file_reader as polyhedron_file_reader

import math

class ExplicitStrategy(BaseExplicitStrategy):

    def __init__(self, all_model_parts, creator_destructor, dem_fem_search, DEM_parameters, procedures):

        BaseExplicitStrategy.__init__(self, all_model_parts, creator_destructor, dem_fem_search, DEM_parameters, procedures)
        
        #self.DEM_parameters = DEM_parameters
        #if "PostSkinSphere" in DEM_parameters.keys():
        #    self.print_skin_sphere = DEM_parameters["PostSkinSphere"].GetBool()
        self.polyhedron_model_part = all_model_parts.Get("PolyhedronPart")


    def SetVariablesAndOptions(self):

        BaseExplicitStrategy.SetVariablesAndOptions(self)

        if "PolyhedronFileName" in self.DEM_parameters.keys():
            polyhedron_file_name = self.DEM_parameters["PolyhedronFileName"].GetString()
            [name, list_of_vertices_list, list_of_faces_list, list_of_size, list_of_volume] = polyhedron_file_reader.ReadPolyhedronFile(polyhedron_file_name)
            pre_utils = PreUtilities()
            for properties in self.spheres_model_part.Properties:
                pre_utils.SetPolyhedronInformationInProperties(name, list_of_vertices_list, list_of_faces_list, list_of_size, list_of_volume, properties)
        
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

    def AddDofs(self):
        # this can safely be called also for restarts, it is internally checked if the dofs exist already
        spheres_model_part = self.all_model_parts.Get("SpheresPart")
        dem_inlet_model_part = self.all_model_parts.Get("DEMInletPart")
        cluster_model_part = self.all_model_parts.Get("ClusterPart")
        polyhedron_model_part = self.all_model_parts.Get("PolyhedronPart")

        model_part_list = [spheres_model_part, cluster_model_part, dem_inlet_model_part, polyhedron_model_part]
        variable_list = [VELOCITY_X, VELOCITY_Y, VELOCITY_Z, ANGULAR_VELOCITY_X, ANGULAR_VELOCITY_Y, ANGULAR_VELOCITY_Z]
        for model_part in model_part_list:
            for variable in variable_list:
                VariableUtils().AddDof(variable, model_part)
            self.Procedures.KratosPrintInfo("DOFs for the DEM solution added correctly")

    def InitializeSolutionStep(self):
        time = self.spheres_model_part.ProcessInfo[TIME]
        self.FixDOFsManually(time)
        self.cplusplus_strategy.ResetPrescribedMotionFlagsRespectingImposedDofs()
        self.cplusplus_strategy.ResetPrescribedMotionFlagsRespectingImposedDofsForPolyhedron()
        self.FixExternalForcesManually(time)

        self.cplusplus_strategy.InitializeSolutionStep()

    def RebuildListOfPolyhedronParticles(self):
        self.cplusplus_strategy.RebuildListOfPolyhedronParticles()

    def SetNormalRadiiOnAllParticles(self):
        #self.cplusplus_strategy.SetNormalRadiiOnAllParticles(self.polyhedron_model_part)
        pass

    def AdvanceInTime(self, time):
        """This function updates and return the current simulation time
        """
        time += self.dt
        self._UpdateTimeInModelParts(time)

        return time

    def _UpdateTimeInModelParts(self, time, is_time_to_print = False):
        spheres_model_part = self.all_model_parts.Get("SpheresPart")
        cluster_model_part = self.all_model_parts.Get("ClusterPart")
        dem_inlet_model_part = self.all_model_parts.Get("DEMInletPart")
        rigid_face_model_part = self.all_model_parts.Get("RigidFacePart")
        polyhedron_model_part = self.all_model_parts.Get("PolyhedronPart")

        self._UpdateTimeInOneModelPart(spheres_model_part, time, self.dt, is_time_to_print)
        self._UpdateTimeInOneModelPart(cluster_model_part, time, self.dt, is_time_to_print)
        self._UpdateTimeInOneModelPart(dem_inlet_model_part, time, self.dt, is_time_to_print)
        self._UpdateTimeInOneModelPart(rigid_face_model_part, time, self.dt, is_time_to_print)
        self._UpdateTimeInOneModelPart(polyhedron_model_part, time, self.dt, is_time_to_print)

    @classmethod
    def _UpdateTimeInOneModelPart(self, model_part, time, dt, is_time_to_print = False):
        ''' This method is redirected to its cpp version with improved speed.
        It also has been updated to classmethod and args so it can be called from external App
        '''

        AuxiliaryUtilities().UpdateTimeInOneModelPart(model_part, time, dt, is_time_to_print)
