##############################################################
# Author: Andres Florez, Universidad de los Andes, Colombia. #
##############################################################

ObjSuf = o
SrcSuf = cc
ExeSuf =
DllSuf = so
OutPutOpt = -o
HeadSuf = h

ROOTCFLAGS = $(shell root-config --cflags) -O2
ROOTLIBS = $(shell root-config --libs) -O2

# Linux with egcs

CXX = g++ -std=c++11
CXXFLAGS += $(ROOTCFLAGS) -I./src/ -Wall 

LD = g++ -std=c++11
LDFLAGS += $(ROOTLIBS) -Wall -l TreePlayer -lboost_program_options -Wno-unused-but-set-variable

SOFLAGS = -shared
LIBS =

SRCDIR = src
OBJDIR = obj

#------------------------------------------------------------------------------
SOURCES = $(wildcard $(SRCDIR)/*.cc)
OBJECTS = $(SOURCES:$(SRCDIR)/%.cc=$(OBJDIR)/%.o)
#------------------------------------------------------------------------------

all: Analyzer_tautau Analyzer_tautau_data

Analyzer_tautau: $(OBJECTS) obj/Analyzer_tautau.O
	$(LD) -o $@ $^ $(LDFLAGS) $(LIBS)	

$(OBJDIR)/%.o: $(SRCDIR)/%.cc
	$(CXX) -o $@  -c $< $(CXXFLAGS)

$(OBJDIR)/%.O: $(SRCDIR)/%.C
	$(CXX) -o $@  -c $< $(CXXFLAGS)


%: $(OBJDIR)/%.o $(OBJDIR)/%.O
	$(LD) -o $@ $< $(LIBS) $(LDFLAGS)

clean:
	@echo "Cleaning..."
	@ls $(OBJDIR)
	@rm -f $(OBJECTS) obj/*.O

.SUFFIXES: .$(SrcSuf) .cc .o .so

Analyzer_tautau_data: $(OBJECTS) obj/Analyzer_tautau_data.O
	$(LD) -o $@ $^ $(LDFLAGS) $(LIBS)

$(OBJDIR)/%.o: $(SRCDIR)/%.cc
	$(CXX) -o $@  -c $< $(CXXFLAGS)

$(OBJDIR)/%.O: $(SRCDIR)/%.C
	$(CXX) -o $@  -c $< $(CXXFLAGS)


%: $(OBJDIR)/%.o $(OBJDIR)/%.O
	$(LD) -o $@ $< $(LIBS) $(LDFLAGS)

clean:
	@echo "Cleaning..."
	@ls $(OBJDIR)
	@rm -f $(OBJECTS) obj/*.O

.SUFFIXES: .$(SrcSuf) .cc .o .so
