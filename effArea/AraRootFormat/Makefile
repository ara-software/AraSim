#############################################################################
## Makefile -- New Version of my Makefile that works on both linux
##              and mac os x
## Ryan Nichol <rjn@hep.ucl.ac.uk>
##############################################################################
include ../StandardDefinitions.mk


#Generic and Site Specific Flags
CXXFLAGS     += $(ROOTCFLAGS)  $(SYSINCLUDES) -I. #$(FFTFLAG)
LDFLAGS      += -g $(ROOTLDFLAGS) 

LIBS          = $(ROOTLIBS) -lMinuit $(SYSLIBS) #$(FFTLIBS) -lMathMore 
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)

#Now the bits we're actually compiling
ROOT_LIBRARY = libAraSimEvent.a
LIB_OBJS =  UsefulAraStationEvent.o UsefulIcrrStationEvent.o #araRootDict.o
CLASS_HEADERS = UsefulAraStationEvent.h UsefulIcrrStationEvent.h araIcrrDefines.h


all : $(ROOT_LIBRARY) Makefile.arch

Makefile.arch:
	@cp ../Makefile.arch ./ $

#The library
$(ROOT_LIBRARY) : $(LIB_OBJS) 
	@echo "Linking $@ ..."
ifeq ($(PLATFORM),macosx)
# We need to make both the .dylib and the .so
	 	$(LD) $(SOFLAGS)$@ $(LDFLAGS) $(LIBS) $^ $(OutPutOpt) $@
 ifneq ($(subst $(MACOSX_MINOR),,123456),123456)
 ifeq ($(MACOSX_MINOR),4)
 		$(LD) -sf $@ $(subst .$(DllSuf),.so,$@)
else
		$(LD) $(LDFLAGS) $^ -bundle -undefined $(UNDEFOPT) $(OutPutOpt) $(subst .$(DllSuf),.so,$@)
 endif
 endif
 else
	$(LD) $(SOFLAGS) $(LDFLAGS) $(LIBS) $(LIB_OBJS) -o $@
 endif


%.o: %.cxx
	@echo "<**Compiling**> "$<
	$(CXX) $(CXXFLAGS) -c $< -o  $@

%.o: %.C
	@echo "<**Compiling**> "$<
	$(CXX) $(CXXFLAGS) -c $< -o  $@
 
#%.$(OBJSUF) : %.$(SRCSUF)
#	@echo "<**Compiling**> "$<
#	$(CXX) $(CXXFLAGS) -c $< -o  $@

araRootDict.C: $(CLASS_HEADERS)
	@echo "Generating dictionary ..."
	@ rm -f *Dict* 
	rootcint $@ -c $(CLASS_HEADERS) LinkDef.h

install:
	install -d ../include
	install -d ../lib
	install -c -m 755 $(ROOT_LIBRARY) ../lib
	install -c -m 644 $(CLASS_HEADERS) ../include 

clean:
	@rm -f *Dict*
	@rm -f *.${OBJSUF}
	@rm -f $(LIBRARY)
	@rm -f $(ROOT_LIBRARY)
	@rm -f $(subst .$(DLLSUF),.so,$(ROOT_LIBRARY))	
	@rm -f $(TEST)
	@rm -f ../lib/$(ROOT_LIBRARY)
#############################################################################



