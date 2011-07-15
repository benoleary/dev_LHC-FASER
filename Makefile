#
# Examples Makefile.
#

SHELL = /bin/sh



# There is no default behaviour, so remind user.
all:
	$(CXX) $(CXXFLAGS) dev_LHC-FASER_unit_tester_main.cpp LHC-FASER_base_lepton_distribution_stuff.cpp LHC-FASER_cross-section_stuff.cpp LHC-FASER_derived_lepton_distributions.cpp LHC-FASER_global_stuff.cpp LHC-FASER_input_handling_stuff.cpp LHC-FASER_sparticle_decay_stuff.cpp ./CppSLHA/*.cpp -o LHC-FASER_unit_tester_main.exe &> compileoutput.txt\
	
