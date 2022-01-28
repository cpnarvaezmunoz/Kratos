from __future__ import print_function, absolute_import, division

import KratosMultiphysics #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

def Wait():
	input("Press Something")

#including kratos path
import sys
import math
import time
from KratosMultiphysics import *
from KratosMultiphysics.EhdtestApplication import *
#import KratosMultiphysics    #we import the KRATOS  
import KratosMultiphysics.EhdtestApplication as Poisson       #and now our application. note that we can import as many as we need to solve our specific problem 
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities


import KratosMultiphysics.EhdtestApplication.poisson_solver_ehd as poisson_solver_ehd           #we import the python file that includes the commands that we need
#import KratosMultiphysics.EhdtestApplication.ElectricfieldCalculateProcess as ElectricfieldCalculateProcess
#setting the domain size for the problem to be solved
domain_size = 2  # 2D problem  

#defining a model part
model = Model()
model_part = model.CreateModelPart("Main");  #we create a model part  
###################################################################################
model_part.AddNodalSolutionStepVariable(Poisson.EPOTENCIAL)
#model_part.Properties[0].SetValue(Poisson.PERMITTIVITY)
model_part.AddNodalSolutionStepVariable(Poisson.PCHARGE)
model_part.AddNodalSolutionStepVariable(Poisson.SCHARGE)
model_part.AddNodalSolutionStepVariable(Poisson.VCHARGE)
model_part.AddNodalSolutionStepVariable(Poisson.EDISPLACEMENT)
model_part.AddNodalSolutionStepVariable(Poisson.EFIELD)   #electric field
model_part.AddNodalSolutionStepVariable(Poisson.ICOEFFICIENT)
model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
#####################################################################################
# Function to define the initial water body by setting a distance field
# Negative distance value represent water
# Positive values represent the air-filled volume

parameter_file = open("ProjectParameters.json",'r')  # we read the ProjectParameters.json file
custom_settings = Parameters(parameter_file.read())

poisson_solver_ehd = poisson_solver_ehd.CreateSolver(model, custom_settings)
poisson_solver_ehd.AddVariables()  #from the static_poisson_solver.py we call the function Addvariables so that the model part we have just created has the needed variables 

#
 # (note that our model part does not have nodes or elements yet) 

 #now we proceed to use the GID interface (both to import the infomation inside the .mdpa file and later print the results in a file  
gid_mode = GiDPostMode.GiD_PostAscii  #we import the python file that includes the commands that we need  
multifile = MultiFileFlag.SingleFile
deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
write_conditions = WriteConditionsFlag.WriteElementsOnly
gid_io = GidIO("art4", gid_mode,multifile, deformed_mesh_flag, write_conditions)

model_part_io = ModelPartIO("example")          # we set the name of the .mdpa file  
model_part_io.ReadModelPart(model_part)         # we load the info from the .mdpa 

# we create a mesh for the postprocess  
mesh_name = 0.0
gid_io.InitializeMesh( mesh_name )
gid_io.WriteMesh((model_part).GetMesh())
gid_io.FinalizeMesh()

#the buffer size should be set up here after the mesh is read for the first time  (this is important for transcient problems, in this static problem =1 is enough)  
model_part.SetBufferSize(1)

# we add the DoFs  
poisson_solver_ehd.AddDofs()

# Initialize the solver
#electricforce_calculate_process = Poisson.ElectricforceCalculateProcess(model_part)
poisson_solver_ehd.Initialize()

#def ModifyInitialGeometry(self):
Yc = 0.001	# plene cut 
#Xc = 5.0e-1   		# drop
#Yc = 5.0e-1		# 

for node in model_part.Nodes:   
        dy = Yc + node.Y
        distance = dy
        #dx = (node.X - Xc)
        #dy = (node.Y - Yc)
        #distance = (math.sqrt( (dx)**2 + (dy)**2) - 0.25**2)
      
        node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, distance)
        if node.Id == 4556:
            node.SetSolutionStepValue(Poisson.EPOTENCIAL,0,0.0)
            node.Fix(Poisson.EPOTENCIAL)

print ("about to solve!")    
poisson_solver_ehd.Solve()
#electricforce_calculate_process.Execute()

#KratosMultiphysics.ComputeNonHistoricalNodalGradientProcess(
#    model_part,
#    ep).Execute()
#print (ep)
print ("Solved!")  

#and we print the results  
gid_io.InitializeResults(mesh_name,(model_part).GetMesh()) 
gid_io.WriteNodalResults(Poisson.EPOTENCIAL,model_part.Nodes,0,0)
#gid_io.PrintOnGaussPoints(Poisson.EFORCEPOS,model_part,0)
#gid_io.PrintOnGaussPoints(Poisson.EFORCENEG,model_part,0)
gid_io.PrintOnGaussPoints(Poisson.EFIELD,model_part,0)
gid_io.WriteNodalResults(KratosMultiphysics.DISTANCE,model_part.Nodes,0,0)
gid_io.FinalizeResults()

# since we have already calculated the temp, we can get the mean value
# first the constructor (it could have been called before)
#try:
#	calc_mean = Poisson.ElectricfieldCalculateProcess(model_part)
    # and we calculate!
#	calc_mean.Execute()
#except:
#	pass
