import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW
from KratosMultiphysics.ShallowWaterApplication.utilities.wave_factory import WaveTheoryFactory


def Factory(settings, Model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return WaveGeneratorProcess(Model, settings["Parameters"])


class WaveGeneratorProcess(KM.Process):

    __formulation = {
        "primitive_variables" : SW.Formulation.PrimitiveVariables,
        "conserved_variables" : SW.Formulation.ConservativeVariables
    }


    @staticmethod
    def GetDefaultParameters():
        """Default parameters for wave generator process.

        A zero value in the wave specifications will be ignored.
        The direction can be specified with a vector or the 'normal' keyword.
        """
        return KM.Parameters("""{
            "model_part_name"          : "model_part",
            "formulation"              : "primitive_variables",
            "interval"                 : [0.0, "End"],
            "direction"                : [0.0, 0.0, 0.0],
            "normal_positive_outwards" : true,
            "smooth_time"              : 0.0,
            "wave_specifications"      : {
                "wave_theory"               : "boussinesq",
                "period"                    : 0.0,
                "wavelength"                : 0.0,
                "amplitude"                 : 0.0,
                "get_depth_from_model_part" : true
            }
        }""")


    def __init__(self, model, settings ):
        """"""
        KM.Process.__init__(self)

        # Check if the direction is specified by 'normal'
        self.settings = settings
        self.direction_by_normal = False
        if self.settings.Has("direction"):
            if self.settings["direction"].IsString():
                if self.settings["direction"].GetString() == "normal":
                    self.settings["direction"].SetVector(KM.Vector(3))
                    self.direction_by_normal = True
        else:
            self.direction_by_normal = True

        # Overwrite the default settings with user-provided parameters
        self.settings.ValidateAndAssignDefaults(self.GetDefaultParameters())

        # Get the custom settings
        self.model_part = model[self.settings["model_part_name"].GetString()]
        self.interval = KM.IntervalUtility(self.settings)
        self.formulation = self.__formulation[self.settings["formulation"].GetString()]
        self.fix_dofs = True


    def ExecuteInitialize(self):
        """Initialize the wave theory and the periodic function."""

        # Check the direction
        if self.direction_by_normal:
            direction = self._CalculateUnitNormal()
            if self.settings["normal_positive_outwards"].GetBool():
                direction *= -1
            self.settings["direction"].SetVector(direction)

        # Setup the wave theory
        wave_settings = self.settings["wave_specifications"]
        wave_theory = WaveTheoryFactory(self.model_part, wave_settings)

        # Creation of the parameters for the c++ process
        velocity_parameters = KM.Parameters("""{}""")
        velocity_parameters.AddEmptyValue("amplitude").SetDouble(wave_theory.horizontal_velocity)
        velocity_parameters.AddEmptyValue("wavelength").SetDouble(wave_theory.wavelength)
        velocity_parameters.AddEmptyValue("period").SetDouble(wave_theory.period)
        velocity_parameters.AddValue("phase", self.settings["wave_specifications"]["t_shift"])
        velocity_parameters.AddValue("shift", self.settings["wave_specifications"]["x_shift"])
        velocity_parameters.AddValue("smooth_time", self.settings["smooth_time"])
        velocity_parameters.AddValue("direction", self.settings["direction"])

        self.velocity_process = SW.ApplySinusoidalFunctionToVector(self.model_part, KM.VELOCITY, velocity_parameters)


    def Check(self):
        self.velocity_process.Check()


    def ExecuteInitializeSolutionStep(self):
        if self._IsInInterval():
            self.velocity_process.ExecuteInitializeSolutionStep()
            if self.formulation == SW.Formulation.ConservativeVariables:
                SW.ShallowWaterUtilities().ComputeMomentum(self.model_part)
            if self.fix_dofs:
                if self.formulation == SW.Formulation.PrimitiveVariables:
                    KM.VariableUtils().ApplyFixity(KM.VELOCITY_X, True, self.model_part.Nodes)
                    KM.VariableUtils().ApplyFixity(KM.VELOCITY_Y, True, self.model_part.Nodes)
                if self.formulation == SW.Formulation.ConservativeVariables:
                    KM.VariableUtils().ApplyFixity(KM.MOMENTUM_X, True, self.model_part.Nodes)
                    KM.VariableUtils().ApplyFixity(KM.MOMENTUM_Y, True, self.model_part.Nodes)


    def _IsInInterval(self):
        """Returns if we are inside the time interval or not."""
        current_time = self.model_part.ProcessInfo[KM.TIME]
        return self.interval.IsInInterval(current_time)


    def _CalculateUnitNormal(self):
        geometry = self.model_part.Conditions.__iter__().__next__().GetGeometry()
        normal = geometry.UnitNormal()
        return normal
