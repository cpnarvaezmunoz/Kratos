import KratosMultiphysics
import KratosMultiphysics.ShallowWaterApplication

def run():
    KratosMultiphysics.Tester.SetVerbosity(KratosMultiphysics.Tester.Verbosity.FAILED_TESTS_OUTPUTS)
    KratosMultiphysics.Tester.RunTestSuite("ShallowWaterApplicationFastSuite")

if __name__ == '__main__':
    run()
