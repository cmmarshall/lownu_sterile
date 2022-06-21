
CXX = g++
ROOTFLAGS = `root-config --libs --cflags`

# ROOT Flags are incomplete.
#LDLIBS += -L$(ROOTSYS)/lib
#INCLUDE += -I$(ROOTSYS)/include

all: DoTemplateFit

DoTemplateFit : DoTemplateFit.cxx
	$(CXX) -o DoTemplateFit DoTemplateFit.cxx $(ROOTFLAGS) 

clean:
	-rm -f DoTemplateFit
	-rm -f DoTemplateFit.o

remake:
	make clean
	make all
