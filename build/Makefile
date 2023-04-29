ROOT 		= ..
CXX  		= mpicxx
#CXX  		= cxx
#HCXX  	= h5cc
CXXFLAGS = -g -I$(ROOT)/include #-std=c++11
#LIBS		= -L$(ROOT)/lib -lm -L/usr/lib64/
#LIBS		= -L$(ROOT)/lib -lm -L/usr/lib64/ -lnetcdf_c++
#HLIBS		= -L/home/user/xfer-software/hdf5/hdf5-1.10.5/lib -lhdf5
LIBS		= -lhdf5 -L$(ROOT)/lib -lm -lz -lsz -lmpi -ldl -L/usr/lib64 #$(HLIBS)
OBJS		= read_cmdline.o readHDF5c.o writeHDF5c.o check_for_failure.o spf_communication.o voxel_dynamics.o read_parameter_file.o rand1D.o jumps.o rungekutta.o stochastic_rates.o free_energy_gl.o #ensemble_quantities.o
SRCS		= ../src/read_cmdline.cpp ../src/check_for_failure.cpp ../src/readHDF5c.cpp ../src/writeHDF5c.cpp ../src/spf_communication.cpp ../src/voxel_dynamics.cpp ../src/read_parameter_file.cpp ../src/rand1D.cpp ../src/jumps.cpp ../src/rungekutta.cpp ../src/stochastic_rates.cpp ../src/free_energy_gl.cpp #../src/ensemble_quantities.cpp
#HSRCS		= ../src/check_for_failure.cpp ../src/readHDF5c.cpp ../src/writeHDF5c.cpp
TARGET	= $(ROOT)/bin/test.x

#random:		$(OBJS) $(HOBJS)
#	$(CXX) $(ROOT)/src/test_random_processes.cpp $(OBJS) $(LIBS) -o $(ROOT)/bin/test_random_processes.x
#	#$(CXX) $(ROOT)/src/test_random_processes.cpp $(OBJS) $(HOBJS) $(LIBS) -o $(TARGET)
#
#poisson: $(OBJS) $(HOBJS)
#	$(CXX) $(ROOT)/src/test_poisson_process.cpp $(OBJS) $(HOBJS) $(LIBS) -o $(ROOT)/bin/test_poisson.x
#
marcus: $(OBJS) $(HOBJS)
	$(CXX) $(ROOT)/src/test_marcus.cpp $(OBJS) $(HOBJS) $(LIBS) -o $(ROOT)/bin/test_marcus.x
#
drift: $(OBJS) $(HOBJS)
	$(CXX) $(ROOT)/src/test_marcus_drift.cpp $(OBJS) $(HOBJS) $(LIBS) -o $(ROOT)/bin/test_marcus_drift.x

spf: $(OBJS) #$(HOBJS)
	$(CXX) $(ROOT)/src/spf_pure_gaussian_single.cpp $(OBJS) $(CXXFLAGS) $(LIBS) -o $(ROOT)/bin/spf_pure_gaussian_single.x
	$(CXX) $(ROOT)/src/spf_pure_gaussian_single_milstein.cpp $(OBJS) $(CXXFLAGS) $(LIBS) -o $(ROOT)/bin/spf_pure_gaussian_single_milstein.x
	$(CXX) $(ROOT)/src/spf_pure_gaussian_separate.cpp $(OBJS) $(CXXFLAGS) $(LIBS) -o $(ROOT)/bin/spf_pure_gaussian_separate.x
	$(CXX) $(ROOT)/src/spf_pure_gaussian_separate_milstein.cpp $(OBJS) $(CXXFLAGS) $(LIBS) -o $(ROOT)/bin/spf_pure_gaussian_separate_milstein.x
	$(CXX) $(ROOT)/src/spf_pure_poisson_simple.cpp $(OBJS) $(CXXFLAGS) $(LIBS) -o $(ROOT)/bin/spf_pure_poisson_simple.x
	$(CXX) $(ROOT)/src/spf_pure_gaussian_simple.cpp $(OBJS) $(CXXFLAGS) $(LIBS) -o $(ROOT)/bin/spf_pure_gaussian_simple.x

spf-poisson-simple: $(OBJS)
	$(CXX) $(ROOT)/src/spf_pure_poisson_simple.cpp $(OBJS) $(CXXFLAGS) $(LIBS) -o $(ROOT)/bin/spf_pure_poisson_simple.x

spf-gaussian-simple: $(OBJS)
	$(CXX) $(ROOT)/src/spf_pure_gaussian_simple.cpp $(OBJS) $(CXXFLAGS) $(LIBS) -o $(ROOT)/bin/spf_pure_gaussian_simple.x

spf-gaussian-pairwise-simple: $(OBJS)
	$(CXX) $(ROOT)/src/spf_pure_gaussian_pairwise_simple.cpp $(OBJS) $(CXXFLAGS) $(LIBS) -o $(ROOT)/bin/spf_pure_gaussian_pairwise_simple.x

spf-poisson-pairwise-simple: $(OBJS)
	$(CXX) $(ROOT)/src/spf_pure_poisson_pairwise_simple.cpp $(OBJS) $(CXXFLAGS) $(LIBS) -o $(ROOT)/bin/spf_pure_poisson_pairwise_simple.x

spf-poisson-pairwise-drift-simple: $(OBJS)
	$(CXX) $(ROOT)/src/spf_pure_poisson_pairwise_drift_simple.cpp $(OBJS) $(CXXFLAGS) $(LIBS) -o $(ROOT)/bin/spf_pure_poisson_pairwise_drift_simple.x

spf-gaussian-pairwise-double-well: $(OBJS)
	$(CXX) $(ROOT)/src/spf_pure_gaussian_pairwise_double_well.cpp $(OBJS) $(CXXFLAGS) $(LIBS) -o $(ROOT)/bin/spf_pure_gaussian_pairwise_double_well.x

spf-poisson-pairwise-double-well: $(OBJS)
	$(CXX) $(ROOT)/src/spf_pure_poisson_pairwise_double_well.cpp $(OBJS) $(CXXFLAGS) $(LIBS) -o $(ROOT)/bin/spf_pure_poisson_pairwise_double_well.x

spf-poisson-pairwise-drift-double-well: $(OBJS)
	$(CXX) $(ROOT)/src/spf_pure_poisson_pairwise_drift_double_well.cpp $(OBJS) $(CXXFLAGS) $(LIBS) -o $(ROOT)/bin/spf_pure_poisson_pairwise_drift_double_well.x

spf-poisson-double-well: $(OBJS)
	$(CXX) $(ROOT)/src/spf_pure_poisson_double_well.cpp $(OBJS) $(CXXFLAGS) $(LIBS) -o $(ROOT)/bin/spf_pure_poisson_double_well.x

spf-gaussian-double-well: $(OBJS)
	$(CXX) $(ROOT)/src/spf_pure_gaussian_double_well.cpp $(OBJS) $(CXXFLAGS) $(LIBS) -o $(ROOT)/bin/spf_pure_gaussian_double_well.x

spf-gaussian-double-well-ito: $(OBJS)
	$(CXX) $(ROOT)/src/spf_pure_gaussian_double_well_ito.cpp $(OBJS) $(CXXFLAGS) $(LIBS) -o $(ROOT)/bin/spf_pure_gaussian_double_well_ito.x

spf-srscd-3d: $(OBJS)
	$(CXX) $(ROOT)/src/spf_srscd_3d.cpp $(OBJS) $(CXXFLAGS) $(LIBS) -o $(ROOT)/bin/spf_srscd_3d.x


gl: $(OBJS)
	$(CXX) $(ROOT)/src/spf_ginzberg_landau_3d.cpp $(OBJS) $(CXXFLAGS) $(LIBS) -o $(ROOT)/bin/spf_ginzberg_landau_3d.x

gl-deterministic: $(OBJS)
	$(CXX) $(ROOT)/src/spf_ginzberg_landau_3d_deterministic.cpp $(OBJS) $(CXXFLAGS) $(LIBS) -o $(ROOT)/bin/spf_ginzberg_landau_3d_deterministic.x

gl-marcus: $(OBJS)
	$(CXX) $(ROOT)/src/spf_ginzberg_landau_3d_marcus.cpp $(OBJS) $(CXXFLAGS) $(LIBS) -o $(ROOT)/bin/spf_ginzberg_landau_3d_marcus.x

stats: $(OBJS) #$(HOBJS)
	$(CXX) $(ROOT)/src/statsHDF5.cpp $(OBJS) $(CXXFLAGS) $(LIBS) -o $(ROOT)/bin/spf_stats.x

1d-walks: $(OBJS)
	$(CXX) $(ROOT)/src/walk_ensembles_1d.cpp  $(OBJS) $(CXXFLAGS) $(LIBS) -o $(ROOT)/bin/spf_walk_ensembles_1d.x

## record Poisson and Brownian random walks in 1-D with average of 0
1d-walks-mu0: $(OBJS)
	$(CXX) $(ROOT)/src/walk_ensembles_mu0_1d.cpp  $(OBJS) $(CXXFLAGS) $(LIBS) -o $(ROOT)/bin/spf_walk_ensembles_mu0_1d.x

1d-escape: $(OBJS)
	$(CXX) $(ROOT)/src/escape_time_ensembles_1d.cpp  $(OBJS) $(CXXFLAGS) $(LIBS) -o $(ROOT)/bin/spf_escape_time_ensembles_1d.x

1d-exponential: $(OBJS)
	$(CXX) $(ROOT)/src/spf_convergencetoexponential_1d.cpp  $(OBJS) $(CXXFLAGS) $(LIBS) -o $(ROOT)/bin/spf_convergencetoexponential_1d.x

#$(CXX) $(ROOT)/src/spf_pure_gaussian.cpp $(OBJS) $(HOBJS) $(LIBS) -o $(ROOT)/bin/spf_pure_gaussian.x

$(OBJS):
	$(CXX) -g -c $(CXXFLAGS) $(LIBS) $(SRCS)

#$(HOBJS):
#	$(CXX) -g -c $(CXXFLAGS) $(HLIBS) $(HSRCS)

all:	$(TARGET)
#all:	poisson random marcus

clean:
	rm -f $(OBJS) $(HOBJS) $(TARGET)  ../bin/spf_1d.x
