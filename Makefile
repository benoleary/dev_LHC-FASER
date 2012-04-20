#
# Example Makefile.
#

SHELL = /bin/sh

CC=g++
CFLAGS=-c -Wall
LDFLAGS=
SOURCES=LHC-FASER.cpp \
 LHC-FASER_base_electroweak_cascade_stuff.cpp \
 LHC-FASER_base_kinematics_stuff.cpp \
 LHC-FASER_base_lepton_distribution_stuff.cpp \
 LHC-FASER_charged_electroweak_cascade_stuff.cpp \
 LHC-FASER_cross-section_stuff.cpp \
 LHC-FASER_derived_lepton_distributions.cpp \
 LHC-FASER_electroweak_cascade_collection_stuff.cpp \
 LHC-FASER_full_cascade_stuff.cpp \
 LHC-FASER_global_stuff.cpp \
 LHC-FASER_input_handling_stuff.cpp \
 LHC-FASER_jet_kinematics_stuff.cpp \
 LHC-FASER_lepton_kinematics_stuff.cpp \
 LHC-FASER_neutral_electroweak_cascade_stuff.cpp \
 LHC-FASER_signal_calculator_stuff.cpp \
 LHC-FASER_signal_data_collection_stuff.cpp \
 LHC-FASER_sparticle_decay_stuff.cpp \
 CppSLHA/CppSLHA_BLOCK.cpp \
 CppSLHA/CppSLHA_EW_scale_spectrum.cpp \
 CppSLHA/CppSLHA_file_reading_stuff.cpp \
 CppSLHA/CppSLHA_global_stuff.cpp \
 CppSLHA/CppSLHA_particle_data.cpp \
 CppSLHA/CppSLHA_PDG_codes_and_data.cpp \
 CppSLHA/CppSLHA_waiting_on_subprocess_executor.cpp \
 CppSLHA/CppSLHA.cpp
OBJECTS=$(SOURCES:.cpp=.o)
TESTEXECUTABLE=dev_LHC-FASER_unit_tester.exe
JONASEXECUTABLE=LHC-FASER_for_JL.exe

all: $(SOURCES) $(TESTEXECUTABLE) $(JONASEXECUTABLE)
	
$(TESTEXECUTABLE): dev_LHC-FASER_unit_tester_main.cpp $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) dev_LHC-FASER_unit_tester_main.cpp -o $@
	
$(JONASEXECUTABLE): LHC-FASER_for_JL.cpp $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) LHC-FASER_for_JL.cpp -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *o $(EXECUTABLE)