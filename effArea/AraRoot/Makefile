#############################################################################
## Makefile -- New Version of my Makefile that works on both linux
##              and mac os x
## Ryan Nichol <rjn@hep.ucl.ac.uk>
##############################################################################
include ../StandardDefinitions.mk


### Customize these two terms to point to system libraries and headers that are not included by default
SYSINCLUDES = 
## e.g -I$(HOME)/tools/RHEL5_x86_64/include
SYSLIBS = 

LD_ARA_UTIL=-L ./
INC_ARA_UTIL=-I ./


#Generic and Site Specific Flags
CXXFLAGS     += $(ROOTCFLAGS)  $(SYSINCLUDES)  -I ./ 
LDFLAGS      += -g $(ROOTLDFLAGS) 


LIBS          = $(ROOTLIBS) -lMathMore -lMinuit $(SYSLIBS) $(LD_ARA_UTIL) -lsqlite3
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)

#Now the bits we're actually compiling
ROOT_LIBRARY = libAraGeom.a
LIB_OBJS = AraAntennaInfo.o AraGeomTool.o AraStationInfo.o eventDictRoot.o

#AraSimTriggerMonitor.o
CLASS_HEADERS = AraAntennaInfo.h AraGeomTool.h  araIcrrStructures.h araIcrrDefines.h AraStationInfo.h 

all : $(ROOT_LIBRARY) Makefile.arch

tests: printStructSizes

Makefile.arch:
	@cp ../Makefile.arch ./ $

#The library
$(ROOT_LIBRARY) : $(LIB_OBJS) 
	@echo "Linking $@ ..."
ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .so
		$(LD) $(SOFLAGS)$@ $(LDFLAGS) $(LIBS) $^ $(OutPutOpt) $@
ifneq ($(subst $(MACOSX_MINOR),,1234),1234)
ifeq ($(MACOSX_MINOR),4)
		ln -sf $@ $(subst .$(DllSuf),.so,$@)
else
		$(LD) -bundle -undefined $(UNDEFOPT) $(LDFLAGS) $^ \
		   $(OutPutOpt) $(subst .$(DllSuf),.so,$@)
endif
endif
else
	$(LD) $(SOFLAGS) $(LDFLAGS) $(LIBS) $(LIB_OBJS) -o $@
endif

%.$(OBJSUF) : %.$(SRCSUF)
	@echo "<**Compiling**> "$<
	$(CXX) $(CXXFLAGS) -c $< -o  $@

%.$(OBJSUF) : %.C
	@echo "<**Compiling**> "$<
	$(CXX) $(CXXFLAGS) $ -c $< -o  $@


eventDictRoot.C: $(CLASS_HEADERS)
	@echo "Generating dictionary ..."
	@ rm -f *Dict* 
	rootcint $@ -c $(CLASS_HEADERS) LinkDef.h

install: $(ROOT_LIBRARY)
ifeq ($(PLATFORM),macosx)
	install -c -m 755 $(ROOT_LIBRARY) $(subst .$(DLLSUF),.so,$(ROOT_LIBRARY)) ../lib
else
	install -c -m 755 $(ROOT_LIBRARY) ../lib
endif		
	install -c -m 644  $(CLASS_HEADERS) ../include


clean:
	@rm -f *Dict*
	@rm -f *.${OBJSUF}
	@rm -f $(LIBRARY)
	@rm -f $(ROOT_LIBRARY)
	@rm -f $(subst .$(DLLSUF),.so,$(ROOT_LIBRARY))	
	@rm -f $(TEST)
#############################################################################



