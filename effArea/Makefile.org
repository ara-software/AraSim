#############################################################################
## Makefile -- New Version of my Makefile that works on both linux
##              and mac os x
## Ryan Nichol <rjn@hep.ucl.ac.uk>
##############################################################################
include Makefile.arch

#Site Specific  Flags
ifeq ($(strip $(BOOST_ROOT)),)
	BOOST_ROOT = /usr/local/include
endif
SYSINCLUDES	= -I/usr/include -I$(BOOST_ROOT)
SYSLIBS         = -L/usr/lib
DLLSUF = ${DllSuf}
OBJSUF = ${ObjSuf}
SRCSUF = ${SrcSuf}

#Now the bits we're actually compiling


#Generic and Site Specific Flags
CXXFLAGS     += $(ROOTCFLAGS) $(SYSINCLUDES)
LDFLAGS      += -g $(ROOTLDFLAGS) 

# copy from ray_solver_makefile (removed -lAra part)
LDFLAGS+=-lboost_program_options -L.

# added for Fortran to C++
G77LDFLAGS+=-lg2c
G77	= g77

LIBS          = $(ROOTLIBS) -lMinuit $(SYSLIBS) 
#LIBS          = $(ROOTLIBS) -lgsl -lMinuit $(SYSLIBS) 
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)

ROOT_LIBRARY = libAra.${DLLSUF}
#LIB_OBJS = AraSim.o Detector.o Event.o Efficiencies.o Trigger.o IceModel.o EarthModel.o eventDict.o
#LIB_OBJS =  Vector.o EarthModel.o IceModel.o Trigger.o Ray.o Tools.o Efficiencies.o Event.o Detector.o Position.o Spectra.o RayTrace.o RayTrace_IceModels.o signal.o eventDict.o Settings.o Primaries.o counting.o
LIB_OBJS =  Vector.o EarthModel.o IceModel.o Trigger.o Ray.o Tools.o Efficiencies.o Event.o Detector.o Position.o Spectra.o RayTrace.o RayTrace_IceModels.o signal.o secondaries.o Settings.o Primaries.o counting.o RaySolver.o convolution_model_Fresnel_xyz.o dgausspkg.o divdifdouble.o vegas.o eventDict.o
CCFILE       =  Vector.cc EarthModel.cc IceModel.cc Trigger.cc Ray.cc Tools.cc Efficiencies.cc Event.cc Detector.cc Spectra.cc Position.cc RayTrace.cc signal.cc secondaries.cc RayTrace_IceModels.cc Settings.cc Primaries.cc counting.cc RaySolver.cc
CLASS_HEADERS = Trigger.h Detector.h Settings.h Spectra.h IceModel.h Primaries.h
#LIB_OBJS = convolution_model_Fresnel_xyz.o dgausspkg.o divdifdouble.o vegas.o

PROGRAMS = araSim

all : $(ROOT_LIBRARY) araSim

araSim : $(ROOT_LIBRARY) AraSim.$(SRCSUF)
	@echo "<**Compiling**> "  
	$(LD) $(CXXFLAGS) $(LDFLAGS) AraSim.$(SRCSUF) $(ROOT_LIBRARY) $(LIBS) -o $@
#--------------------------------------------------
# araSim : $(ROOT_LIBRARY) AraSim.$(SRCSUF)
# 	@echo "<**Compiling**> "  
# 	$(LD) $(CXXFLAGS) $(LDFLAGS) $(LIBS) $(LIB_OBJS) -o $@
#-------------------------------------------------- 

#The library
$(ROOT_LIBRARY) : $(LIB_OBJS) 
	@echo "Linking $@ ..."
ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .so
		$(LD) $(SOFLAGS)$@ $(LDFLAGS) $(G77LDFLAGS) $^ $(OutPutOpt) $@
ifneq ($(subst $(MACOSX_MINOR),,1234),1234)
ifeq ($(MACOSX_MINOR),4)
		ln -sf $@ $(subst .$(DllSuf),.so,$@)
else
		$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $(G77LDFLAGS) $^ \
		   $(OutPutOpt) $(subst .$(DllSuf),.so,$@)
endif
endif
else
	$(LD) $(SOFLAGS) $(LDFLAGS) $(G77LDFLAGS) $(LIBS) $(LIB_OBJS) -o $@
endif

#%.$(OBJSUF) : %.$(SRCSUF)
#	@echo "<**Compiling**> "$<
#	$(CXX) $(CXXFLAGS) -c $< -o  $@



%.$(OBJSUF) : %.C
	@echo "<**Compiling**> "$<
	$(CXX) $(CXXFLAGS) $ -c $< -o  $@

%.$(OBJSUF) : %.cc
	@echo "<**Compiling**> "$<
	$(CXX) $(CXXFLAGS) $ -c $< -o  $@


# added for fortran code compiling
%.$(OBJSUF) : %.f
	@echo "<**Compiling**> "$<
	$(G77) -c $<


eventDict.C: $(CLASS_HEADERS)
	@echo "Generating dictionary ..."
	@ rm -f *Dict* 
	rootcint $@ -c $(CLASS_HEADERS) LinkDef.h

clean:
	@rm -f *Dict*
	@rm -f *.${OBJSUF}
	@rm -f $(LIBRARY)
	@rm -f $(ROOT_LIBRARY)
	@rm -f $(subst .$(DLLSUF),.so,$(ROOT_LIBRARY))	
	@rm -f $(TEST)
#############################################################################



