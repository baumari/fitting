# General Makefile for c, c++
# Sources are must be in ./src
# Headers are must be in ./include
# Other file or directory are made automatically.

INCDIR = $(CURDIR)/include
INCDIR_EXTRA1 = /Users/kento/work/physics/kinema/include
INCDIR_EXTRA2 = 
#OUTLIBDIR = $(CURDIR)/lib
INCLIBDIR = /Users/kento/work/physics/kinema/lib
SRCDIR = $(CURDIR)/src
OBJDIR = $(SRCDIR)/obj
BINDIR = $(CURDIR)/bin
TARGET = $(BINDIR)/fitting
#LIB = $(OUTLIBDIR)/libkinema.so
LIB =
SUFFIX = cxx

CXX = `root-config --cxx`
#CXX = g++
#CXX = clang++

#CC =
######################################################
# You should edit above depending on your environment.
######################################################

ifeq ($(CXX),`root-config --cxx`)
	ROOTLIBS = `root-config --libs`
	CXXFLAGS = `root-config --cflags` \
	-O3 -Wall -Wextra -Wno-unused -Wno-long-long -Wno-unused-command-line-argument \
	-fno-common -I$(INCDIR) -I$(ROOTLIBS) -I$(INCDIR_EXTRA1) -I$(INCDIR_EXTRA2) -fPIC -MMD -MP
	LDFLAGS = `root-config --glibs` -lm -lgsl -lgslcblas -L$(INCLIBDIR) -lkinema
#	LDFLAGS = `root-config --glibs` -lm -lkinema 
else
	CXXFLAGS = -O3 -Wall -Wextra -Wno-unused -Wno-long-long -Wno-unused-command-line-argument \
	-fno-common -I$(INCDIR) -I$(INCDIR_EXTRA1) -I$(INCDIR_EXTRA2) -fPIC -MMD -MP
	LDFLAGS = -lm
#	LDFLAGS = -lm -lkinema
endif

SOURCES = $(wildcard $(SRCDIR)/*.$(SUFFIX))
OBJS = $(addprefix $(OBJDIR)/, $(notdir $(SOURCES:.cxx=.o)))
DEPS = $(OBJS:.o=.d)

all: FORCE $(TARGET)

$(TARGET): $(OBJS) $(LIB)
	-mkdir -p $(BINDIR)
	$(CXX) $(LDFLAGS) -o $@ $^

$(OBJDIR)/%.o: $(SRCDIR)/%.$(SUFFIX)
	-mkdir -p $(OBJDIR)
	$(CXX) $(CXXFLAGS) -o $@ -c $<

$(LIB):
	rm -f $(LIB)
	-mkdir -p $(OUTLIBDIR)
	$(CXX) -shared $(OBJS) $(LDFLAGS) -o $(LIB)

FORCE:
	rm -f $(LIB)
clean:
	-rm -f $(OBJS) $(DEPS) $(TARGET) $(LIB)
	-rm -rf $(BINDIR) $(OBJDIR) $(OUTLIBDIR)

ifneq ($(filter clean,$(MAKECMDGOALS)),clean)
-include $(DEPS)
endif

