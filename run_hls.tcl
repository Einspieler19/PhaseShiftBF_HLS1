
##############################################
# Project settings

# Create a project
open_project	-reset PhaseShiftBeamformer_prj

# The source file and test bench
add_files			PhaseShiftBeamformer.cpp
add_files -tb			PhaseshiftBeamformer_test.cpp
# Specify the top-level function for synthesis
set_top				PhaseShiftBeamformer

###########################
# Solution settings

# Create solution1
open_solution -reset solution1

# Specify a Xilinx device and clock period
# - Do not specify a clock uncertainty (margin)
# - Let the  margin to default to 12.5% of clock period
set_part  {xc7v585t-ffg1157-1}
create_clock -period 10
#set_clock_uncertainty 1.25

# Simulate the C code 
csim_design

# Do not perform any other steps
# - The basic project will be opened in the GUI 
exit

