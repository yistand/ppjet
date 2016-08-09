os = $(shell uname -s)

INCFLAGS      = -I$(ROOTSYS)/include


ifeq ($(os),Linux)
#ly CXXFLAGS      = 
CXXFLAGS      = -fPIC #ly	get complaint from complier and add it 
else
#CXXFLAGS      = -O -fPIC -pipe -Wall -Wno-deprecated-writable-strings -Wno-unused-variable -Wno-unused-private-field -Wno-gnu-static-float-init
## for debugging:
 CXXFLAGS      = -g -O0 -fPIC -pipe -Wall -Wno-deprecated-writable-strings -Wno-unused-variable -Wno-unused-private-field -Wno-gnu-static-float-init
endif

ifeq ($(os),Linux)
LDFLAGS       = -g
LDFLAGSS      = -g --shared 
else
LDFLAGS       = -O -Xlinker -bind_at_load -flat_namespace
LDFLAGSS      = -flat_namespace -undefined suppress
LDFLAGSSS     = -bundle
endif

ifeq ($(os),Linux)
CXX          = g++ 
else
CXX          = clang
endif


ROOTLIBS      = $(shell root-config --libs)

LIBPATH       = $(ROOTLIBS) 


###############################################################################
################### Remake when these headers are touched #####################
###############################################################################

###############################################################################
# standard rules
%.o : %.cxx
	@echo 
	@echo COMPILING
	$(CXX) $(CXXFLAGS) $(INCFLAGS) -c $< -o $@

%  : %.o 
	@echo 
	@echo LINKING
	$(CXX) $(LDFLAGS) $(LIBPATH) $(LIBS) $^ -o $@


###############################################################################
###############################################################################
############################# Main Targets ####################################
###############################################################################
all    : Demo4ClassTofMatch  

Demo4ClassTofMatch : Demo4ClassTofMatch.o ClassTofMatchWeight.o ClassTPCWeight.o




###############################################################################
##################################### MISC ####################################
###############################################################################

clean :
	@echo 
	@echo CLEANING
	rm -vf *.o


