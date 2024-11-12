from KratosMultiphysics import *
from KratosMultiphysics.ConvectionDiffusionSedimentApplication import *

def run():
    Tester.SetVerbosity(Tester.Verbosity.PROGRESS) # TESTS_OUTPUTS
    Tester.RunTestSuite("KratosConvectionDiffusionSedimentFastSuite")

if __name__ == '__main__':
    run()
