ROOT 		= ..
CXX  		= mpicxx
#CXX  		= cxx
#HCXX  	= h5cc
CXXFLAGS = -I$(ROOT)/include #-std=c++11
#LIBS		= -L$(ROOT)/lib -lm -L/usr/lib64/
#LIBS		= -L$(ROOT)/lib -lm -L/usr/lib64/ -lnetcdf_c++
#HLIBS		= -L/home/user/xfer-software/hdf5/hdf5-1.10.5/lib -lhdf5
LIBS		= -lhdf5 -L$(ROOT)/lib -lm -lz -lsz -lmpi -ldl -L/usr/lib64 #$(HLIBS)
OBJS		= read_cmdline.o readHDF5c.o writeHDF5c.o check_for_failure.o spf_communication.o voxel_dynamics.o read_parameter_files.o 
#rand1D.o jumps.o rungekutta.o 
SRCS		= ../src/read_cmdline.cpp ../src/check_for_failure.cpp ../src/readHDF5c.cpp ../src/writeHDF5c.cpp ../src/spf_communication.cpp ../src/voxel_dynamics.cpp ../src/read_parameter_files.cpp
#../src/rand1D.cpp ../src/jumps.cpp ../src/rungekutta.cpp
#HSRCS		= ../src/check_for_failure.cpp ../src/readHDF5c.cpp ../src/writeHDF5c.cpp
TARGET	= $(ROOT)/bin/test.x

#random:		$(OBJS) $(HOBJS)
#	$(CXX) $(ROOT)/src/test_random_processes.cpp $(OBJS) $(LIBS) -o $(ROOT)/bin/test_random_processes.x
#	#$(CXX) $(ROOT)/src/test_random_processes.cpp $(OBJS) $(HOBJS) $(LIBS) -o $(TARGET)
#
#poisson: $(OBJS) $(HOBJS)
#	$(CXX) $(ROOT)/src/test_poisson_process.cpp $(OBJS) $(HOBJS) $(LIBS) -o $(ROOT)/bin/test_poisson.x
#
#marcus: $(OBJS) $(HOBJS)
#	$(CXX) $(ROOT)/src/test_marcus.cpp $(OBJS) $(HOBJS) $(LIBS) -o $(ROOT)/bin/test_marcus.x
#
#drift: $(OBJS) $(HOBJS)
#	$(CXX) $(ROOT)/src/test_marcus_drift.cpp $(OBJS) $(HOBJS) $(LIBS) -o $(ROOT)/bin/test_marcus_drift.x

spf: $(OBJS) #$(HOBJS)
	$(CXX) $(ROOT)/src/spf_pure_gaussian.cpp $(OBJS) $(CXXFLAGS) $(LIBS) -o $(ROOT)/bin/spf_pure_gaussian.x
	$(CXX) $(ROOT)/src/spf_pure_poisson.cpp $(OBJS) $(CXXFLAGS) $(LIBS) -o $(ROOT)/bin/spf_pure_poisson.x

stats: $(OBJS) #$(HOBJS)
	$(CXX) $(ROOT)/src/statsHDF5.cpp $(OBJS) $(CXXFLAGS) $(LIBS) -o $(ROOT)/bin/spf_stats.x

#$(CXX) $(ROOT)/src/spf_pure_gaussian.cpp $(OBJS) $(HOBJS) $(LIBS) -o $(ROOT)/bin/spf_pure_gaussian.x

$(OBJS):
	$(CXX) -g -c $(CXXFLAGS) $(LIBS) $(SRCS)

#$(HOBJS):
#	$(CXX) -g -c $(CXXFLAGS) $(HLIBS) $(HSRCS)

all:	$(TARGET)
#all:	poisson random marcus

clean:
	rm -f $(OBJS) $(HOBJS) $(TARGET) 
