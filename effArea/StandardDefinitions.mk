#############################################################################
## Makefile -- New Version of my Makefile that works on both linux
##              and mac os x
## Ryan Nichol <rjn@hep.ucl.ac.uk>
##############################################################################
#Site Specific  Flags
include Makefile.arch

#If you have things in non standard paths (eg. libRootFftwWrapper) append the appropriate -I or -L lines below
SYSINCLUDES	= 
SYSLIBS         = 

DLLSUF = ${DllSuf}
OBJSUF = ${ObjSuf}
SRCSUF = ${SrcSuf}

#Toggles the FFT functions on and off
#USE_FFT_TOOLS=1

ifdef USE_FFT_TOOLS
FFTLIBS = -lRootFftwWrapper -lfftw3 -lMathMore -lgsl 
FFTFLAG = -DUSE_FFT_TOOLS
else
FFTLIBS =
FFTFLAG =
endif

ifdef USE_GOOGLE_PROFILER
SYSLIBS += -L/home/rjn/thirdParty/lib -lprofiler -ltcmalloc
endif
