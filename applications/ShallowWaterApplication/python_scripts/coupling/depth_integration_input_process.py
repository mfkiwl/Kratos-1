import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW
import KratosMultiphysics.MappingApplication as Mapping
from KratosMultiphysics.kratos_utilities import GenerateVariableListFromInput
from KratosMultiphysics.HDF5Application import import_model_part_from_hdf5_process
from KratosMultiphysics.HDF5Application import single_mesh_temporal_input_process
from os import path, listdir

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return DepthIntegrationInputProcess(model, settings["Parameters"])

class DepthIntegrationInputProcess(KM.OutputProcess):
    """DepthIntegrationInputProcess

    Read the depth integrated values from an HDF5 file and set them as boundary conditions.
    """

    def GetDefaultParameters(self):
        default_parameters = KM.Parameters("""{
            "interface_model_part_name"   : "",
            "input_model_part_name"       : "input_model_part",
            "read_historical_database"    : false,
            "interval"                    : [0.0,"End"],
            "list_of_variables"           : ["MOMENTUM"],
            "list_of_variables_to_fix"    : ["MOMENTUM_X","MOMENTUM_Y"],
            "default_time_after_interval" : null,
            "semi_period_after_interval"  : 1.0,
            "swap_yz_axis"                : false,
            "ignore_vertical_component"   : true,
            "file_settings"               : {}
        }""")
        if self.settings.Has("default_time_after_interval"):
            if self.settings["default_time_after_interval"].IsDouble():
                default_parameters["default_time_after_interval"].SetDouble(0.0)
        return default_parameters

    def __init__(self, model, settings):
        """The constructor of the DepthIntegrationInputProcess."""

        KM.OutputProcess.__init__(self)
        self.settings = settings
        self.settings.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.interface_model_part = model[self.settings["interface_model_part_name"].GetString()]
        self.input_model_part = model.CreateModelPart(self.settings["input_model_part_name"].GetString())
        self.interval = KM.IntervalUtility(self.settings)
        self.variables = GenerateVariableListFromInput(self.settings["list_of_variables"])
        self.variables_to_fix = GenerateVariableListFromInput(self.settings["list_of_variables_to_fix"])

        self.hdf5_import = import_model_part_from_hdf5_process.Factory(self._CreateHDF5Parameters(), model)
        self.hdf5_process = single_mesh_temporal_input_process.Factory(self._CreateHDF5Parameters(), model)
        self._GetInputTimes(self.settings['file_settings'])


    def Check(self):
        '''Check the processes.'''
        self.hdf5_import.Check()
        self.hdf5_process.Check()


    def ExecuteInitialize(self):
        '''Read the input_model_part and set the variables.'''
        self.hdf5_import.ExecuteInitialize()
        self._CheckInputVariables()
        self._CreateMapper()
        self._MapToBoundaryCondition()


    def ExecuteInitializeSolutionStep(self):
        '''Set the variables in the input_model_part at the current time.'''
        current_time = self.interface_model_part.ProcessInfo.GetValue(KM.TIME)
        if self.interval.IsInInterval(current_time):
            self._SetCurrentTime()
            self.hdf5_process.ExecuteInitializeSolutionStep()
            self._CheckInputVariables()
            self._MapToBoundaryCondition()
        else:
            if self.settings["default_time_after_interval"] is not None:
                self._SetDefaultTime()
                self.hdf5_process.ExecuteInitializeSolutionStep()
                self._CheckInputVariables()
                self._MapToBoundaryCondition()
                self._SmoothDefaultValue()


    def _GetInputTimes(self, file_settings):
        # Get all the file names
        file_name = file_settings["file_name"].GetString()
        folder_name = path.dirname(file_name)
        file_names = [f for f in listdir(folder_name) if path.isfile(path.join(folder_name, f))]

        # Find the common parts (prefix and suffix) of the found names and the file pattern
        # The different part is the time, we need to store all the available times
        self.file_pattern = path.basename(file_name)
        self.file_pattern = self.file_pattern.replace('<model_part_name>', self.input_model_part.Name)
        prefix = path.commonprefix([self.file_pattern, file_names[0]])
        suffix = self.file_pattern.replace(prefix, '')
        suffix = path.commonprefix([''.join(reversed(suffix)), ''.join(reversed(file_names[0]))])
        suffix = ''.join(reversed(suffix))
        self.times = []
        for f in file_names:
            f = f.replace(prefix, '')
            f = f.replace(suffix, '')
            self.times.append(float(f))
        self.times.sort()


    def _SetCurrentTime(self):
        current_time = self.interface_model_part.ProcessInfo.GetValue(KM.TIME)
        closest_time = next(filter(lambda x: x>current_time, self.times))
        self.input_model_part.ProcessInfo.SetValue(KM.TIME, closest_time)


    def _SetDefaultTime(self):
        default_time = self.settings["default_time_after_interval"].GetDouble()
        self.input_model_part.ProcessInfo.SetValue(KM.TIME, default_time)


    def _CheckInputVariables(self):
        if self.settings["swap_yz_axis"].GetBool():
            if self.settings["read_historical_database"].GetBool():
                SW.ShallowWaterUtilities().SwapYZComponents(KM.MOMENTUM, self.input_model_part.Nodes)
                SW.ShallowWaterUtilities().SwapYZComponents(KM.VELOCITY, self.input_model_part.Nodes)
            else:
                SW.ShallowWaterUtilities().SwapYZComponentsNonHistorical(KM.MOMENTUM, self.input_model_part.Nodes)
                SW.ShallowWaterUtilities().SwapYZComponentsNonHistorical(KM.VELOCITY, self.input_model_part.Nodes)
        if self.settings["ignore_vertical_component"].GetBool():
            if self.settings["read_historical_database"].GetBool():
                KM.VariableUtils().SetVariableToZero(KM.MOMENTUM_Z, self.input_model_part.Nodes)
                KM.VariableUtils().SetVariableToZero(KM.VELOCITY_Z, self.input_model_part.Nodes)
            else:
                KM.VariableUtils().SetNonHistoricalVariableToZero(KM.MOMENTUM_Z, self.input_model_part.Nodes)
                KM.VariableUtils().SetNonHistoricalVariableToZero(KM.VELOCITY_Z, self.input_model_part.Nodes)


    def _MapToBoundaryCondition(self):
        for variable in self.variables:
            if self.settings["read_historical_database"].GetBool():
                self.mapper.Map(variable, variable)
            else:
                self.mapper.Map(variable, variable, KM.Mapper.FROM_NON_HISTORICAL)

        for variable in self.variables_to_fix:
            KM.VariableUtils().ApplyFixity(variable, True, self.interface_model_part.Nodes)


    def _CreateMapper(self):
        mapper_settings = KM.Parameters("""{
            "mapper_type": "nearest_neighbor",
            "echo_level" : 0,
            "search_settings" : {
                "search_radius" : 0.0
            }
        }""")
        min_point = KM.Point([ 1e6,  1e6,  1e6])
        max_point = KM.Point([-1e6, -1e6, -1e6])
        for node in self.interface_model_part.Nodes:
            for i in range(3):
                point = KM.Point(node)
                min_point[i] = min([min_point[i], point[i]])
                max_point[i] = max([max_point[i], point[i]])
        distance = 1.05 * (max_point - min_point).norm_2()
        mapper_settings["search_settings"]["search_radius"].SetDouble(distance)

        self.mapper = KM.MapperFactory.CreateMapper(
            self.input_model_part,
            self.interface_model_part,
            mapper_settings)


    def _CreateHDF5Parameters(self):
        hdf5_settings = KM.Parameters()
        hdf5_settings.AddValue("model_part_name", self.settings["input_model_part_name"])
        hdf5_settings.AddValue("file_settings", self.settings["file_settings"])
        data_settings = KM.Parameters("""{"list_of_variables" : ["MOMENTUM","VELOCITY","HEIGHT"]}""")
        if self.settings["read_historical_database"].GetBool():
            hdf5_settings.AddValue("nodal_solution_step_data_settings", data_settings)
        else:
            hdf5_settings.AddValue("nodal_data_value_settings", data_settings)
        hdf5_process_settings = KM.Parameters()
        hdf5_process_settings.AddValue("Parameters", hdf5_settings)
        return hdf5_process_settings


    def _SmoothDefaultValue(self):
        elapsed_time = self.interface_model_part.ProcessInfo.GetValue(KM.DELTA_TIME)
        semi_period = self.settings["semi_period_after_interval"].GetDouble()
        for variable in self.variables:
            SW.ShallowWaterUtilities().SmoothHistoricalVariable(
                variable,
                self.interface_model_part.Nodes,
                elapsed_time,
                semi_period)