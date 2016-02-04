os = $(shell uname -s)

#INCFLAGS      = -I$(ROOTSYS)/include -I$(FASTJETDIR)/include -I$(PYTHIA8DIR)/include -I$(STARPICOPATH)
INCFLAGS      = -I$(ROOTSYS)/include -I$(FASTJETDIR)/include -I$(PYTHIA8DIR)/include -I$(PYTHIA8DIR)/include/Pythia8/ -I$(PYTHIA8DIR)/include/Pythia8Plugins/ -I$(STARPICOPATH)
INCFLAGS      += -I./src
INCFLAGS      += -I./include


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

LIBPATH       = $(ROOTLIBS) -L$(FASTJETDIR)/lib -L$(PYTHIA8DIR)/lib -L$(STARPICOPATH)
#ly LIBS          = -lfastjet -lfastjettools -lpythia8  -llhapdfdummy -lTStarJetPico
LIBS          = -lfastjet -lfastjettools -lpythia8  -lTStarJetPico

### Unfolding Test
#INCFLAGS      += -I/home/hep/caines/ly247/Software/RooUnfold-1.1.1/src
#LIBPATH       += -L/home/hep/caines/ly247/Software/RooUnfold-1.1.1
#LIBS          += -lRooUnfold



# for cleanup
SDIR          = src
ODIR          = src/obj
BDIR          = bin


###############################################################################
################### Remake when these headers are touched #####################
###############################################################################
INCS = $(SDIR)/JetAnalyzer.hh $(SDIR)/AjParameters.hh

###############################################################################
# standard rules
$(ODIR)/%.o : $(SDIR)/%.cxx $(INCS)
	@echo 
	@echo COMPILING
	$(CXX) $(CXXFLAGS) $(INCFLAGS) -c $< -o $@

$(BDIR)/%  : $(ODIR)/%.o 
	@echo 
	@echo LINKING
	$(CXX) $(LDFLAGS) $(LIBPATH) $(LIBS) $^ -o $@

###############################################################################

###############################################################################
############################# Main Targets ####################################
###############################################################################
all    : $(BDIR)/PicoJetUnderlyingActivity  \
	lib/libMyJetlib.so 


#$(BDIR)/ppInAuAuAj  \
	 $(BDIR)/RandomCone  \
	 $(BDIR)/TreeWithMc \
	 $(BDIR)/MakeSmallerTrees \
	 lib/libMyJetlib.so \
	 doxy

# lib/libMyJetlib.a \

# $(BDIR)/FollowPythiaAj 
# 	 $(BDIR)/ppInMcAj $(BDIR)/PythiaInMcAj \

$(SDIR)/dict.cxx 		: $(SDIR)/ktTrackEff.hh
	cd $(SDIR); rootcint -f dict.cxx -c -I. ./ktTrackEff.hh

$(ODIR)/dict.o 		: $(SDIR)/dict.cxx
$(ODIR)/ktTrackEff.o 	: $(SDIR)/ktTrackEff.cxx $(SDIR)/ktTrackEff.hh



# $(BDIR)/AreaTest : 	 $(ODIR)/AreaTest.o
# $(BDIR)/AreaTest.o : 	 $(SDIR)/AreaTest.cxx

lib/libMyJetlib.so	: $(ODIR)/JetAnalyzer.o $(ODIR)/dict.o $(ODIR)/ktTrackEff.o
	@echo 
	@echo MAKING LIBRARY
	$(CXX) -shared $(LDFLAGS) $(LIBPATH) $(LIBS) $^ -o $@

#	g++ -flat_namespace -undefined suppress dict.o -L/usr/local/root_v5.32_binary_m64/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lpthread -Wl,-rpath,/usr/local/root_v5.32_binary_m64/lib -lm -ldl -o libKtJet.so
# ar -rcs $@ $^

# lib/libMyJetlib.a	: $(ODIR)/JetAnalyzer.o $(ODIR)/dict.o $(ODIR)/ktTrackEff.o
# 	@echo 
# 	@echo MAKING LIBRARY
# 	ar -rcs $@ $^

$(ODIR)/JetAnalyzer.o 		: $(SDIR)/JetAnalyzer.cxx $(INCS)
$(ODIR)/UnderlyingAna.o 	 	: $(SDIR)/UnderlyingAna.cxx $(INCS) $(SDIR)/UnderlyingAna.hh



#Aj
$(BDIR)/PicoJetUnderlyingActivity		: $(ODIR)/PicoJetUnderlyingActivity.o		$(ODIR)/UnderlyingAna.o	 	lib/libMyJetlib.so
$(BDIR)/STARPythiaJetUnderlyingActivity		: $(ODIR)/STARPythiaJetUnderlyingActivity.o	$(ODIR)/StarPythia.o	$(ODIR)/UnderlyingAna.o	 	lib/libMyJetlib.so
$(BDIR)/ppInAuAuAj 	: $(ODIR)/ppInAuAuAj.o 		$(ODIR)/UnderlyingAna.o	 	lib/libMyJetlib.so
$(BDIR)/ppInMcAj	: $(ODIR)/ppInMcAj.o		$(ODIR)/UnderlyingAna.o	 	lib/libMyJetlib.so
$(BDIR)/PythiaAj	: $(ODIR)/PythiaAj.o 		$(ODIR)/UnderlyingAna.o	 	lib/libMyJetlib.so
$(BDIR)/PythiaInAuAuAj	: $(ODIR)/PythiaInAuAuAj.o 	$(ODIR)/UnderlyingAna.o	 	lib/libMyJetlib.so
$(BDIR)/PythiaInMcAj	: $(ODIR)/PythiaInMcAj.o 	$(ODIR)/UnderlyingAna.o	 	lib/libMyJetlib.so
$(BDIR)/TestCountPythia	: $(ODIR)/TestCountPythia.o 	$(ODIR)/UnderlyingAna.o	 	lib/libMyJetlib.so

$(BDIR)/TreeWithMc      : $(ODIR)/TreeWithMc.o		lib/libMyJetlib.so
$(BDIR)/RandomCone	: $(ODIR)/RandomCone.o		$(ODIR)/UnderlyingAna.o	 	lib/libMyJetlib.so


# helper
$(BDIR)/MakeSmallerTrees	: $(ODIR)/MakeSmallerTrees.o	 	lib/libMyJetlib.so

###############################################################################
##################################### MISC ####################################
###############################################################################


doxy: html/index.html

html/index.html : $(INCS) src/* Doxyfile
#	doxygen
	@echo 
	@echo Updating documentation
	( cat Doxyfile ; echo "QUIET=YES" ) | doxygen -

clean :
	@echo 
	@echo CLEANING
	rm -vf $(ODIR)/*.o
	rm -vf $(BDIR)/*
	rm -vf lib/*


.PHONY : clean doxy
