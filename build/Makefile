#########################################################################
#    SPF - Stochastic Phase Field
#    Copyright (C) 2025 
#    Peng Geng <penggeng@g.ucla.edu>
#    Nicholas Huebner Julian <njulian@ucla.edu>
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License along
#    with this program; if not, write to the Free Software Foundation, Inc.,
#    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
#    See the README file in the top-level SPF directory.
#########################################################################

ROOT 		= ..
SRC_DIR         = $(ROOT)/src
BIN_DIR         = $(ROOT)/bin

CXX  		= mpicxx

CXXFLAGS        = -g -I$(ROOT)/include -I$(SRC_DIR)
LDFLAGS         = -L$(ROOT)/lib -L/usr/lib64 -Wl,-rpath,'$$ORIGIN/../lib'
LDLIBS		= -lhdf5 -lm -lz -lsz -lmpi -ldl

SRCS		= $(SRC_DIR)/check_for_failure.cpp \
                  $(SRC_DIR)/readHDF5c.cpp \
                  $(SRC_DIR)/writeHDF5c.cpp \
                  $(SRC_DIR)/spf_communication.cpp \
                  $(SRC_DIR)/voxel_dynamics.cpp \
                  $(SRC_DIR)/read_parameter_file.cpp \
                  $(SRC_DIR)/stochastic_rates.cpp

OBJS            = $(notdir $(SRCS:.cpp=.o))

TARGET	        = $(BIN_DIR)/spf_3d.x

all: $(TARGET)

$(TARGET): $(OBJS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $(OBJS) $(SRC_DIR)/spf_3d.cpp $(LDFLAGS) $(LDLIBS) -o $@

%.o: $(SRC_DIR)/%.cpp
	$(CXX) -c $(CXXFLAGS) $< -o $@

.PHONY: clean
clean:
	rm -f $(OBJS)
