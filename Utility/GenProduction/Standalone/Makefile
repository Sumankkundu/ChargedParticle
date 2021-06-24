SHELL=/usr/bin/env bash



ROOTCONFIG   := $(ROOTSYS)/bin/root-config
ROOTCFLAGS   := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS  := $(shell $(ROOTCONFIG) --ldflags)
ROOTLIBS     := $(shell $(ROOTCONFIG) --libs)
ROOTGLIBS    := $(shell $(ROOTCONFIG) --glibs)
ROOTVERSION  := $(shell $(ROOTCONFIG) --version)
ROOTMAJORVERSION := $(word 1,$(subst ., ,$(ROOTVERSION)))
ROOTCINT=$(ROOTSYS)/bin/rootcint

CXXFLAGS=-isystem $(shell $(ROOTCONFIG) --incdir) -I$(ROOTSYS)/htmldoc -I. -O2 -g -Wall -Wshadow -W -Woverloaded-virtual -fPIC $(ROOTCFLAGS)
LDFLAGS=$(ROOTLDFLAGS) -L. -Wl,-rpath .
CXX=$(shell $(ROOTCONFIG) --cxx)

OSTYPE=$(shell uname)
ROOTLIBS     := -L$(shell $(ROOTCONFIG) --libdir) -lXMLParser $(ROOTLIBS)


# Install directory prefixes.
#PREFIX_BIN=/home/suman/HEP_Package/pythia8235/bin
#PREFIX_INCLUDE=/home/suman/HEP_Package/pythia8235/include
#PREFIX_LIB=/home/suman/HEP_Package/pythia8235/lib
#PREFIX_SHARE=/home/suman/HEP_Package/pythia8235/share/Pythia8


# PYTHIA configuration file.
# Generated on Wed Nov 13 11:14:36 IST 2019 with the user supplied options:
# --with-hepmc2=/user/local
# --with-root=/home/suman/HEP_Package/ROOT616build
# --with-fastjet3=/user/local
# --with-lhapdf6=/user/local

# Install directory prefixes.
PREFIX_BIN=/home/suman/HEP_Package/pythia8301/bin
PREFIX_INCLUDE=/home/suman/HEP_Package/pythia8301/include
PREFIX_LIB=/home/suman/HEP_Package/pythia8301/lib
PREFIX_SHARE=/home/suman/HEP_Package/pythia8301/share/Pythia8

# Compilation flags (see ./configure --help for further documentation).
CXX=g++
CXX_COMMON=-O2 -std=c++11 -pedantic -W -Wall -Wshadow -fPIC
CXX_SHARED=-shared
CXX_SONAME=-Wl,-soname,
LIB_SUFFIX=.so

EVTGEN_USE=false
EVTGEN_CONFIG=
EVTGEN_BIN=
EVTGEN_INCLUDE=
EVTGEN_LIB=

FASTJET3_USE=true
FASTJET3_CONFIG=fastjet-config
FASTJET3_BIN=/user/local/
FASTJET3_INCLUDE=-I/user/local
FASTJET3_LIB=-L/user/local -Wl,-rpath,/user/local -lfastjet

HEPMC2_USE=true
HEPMC2_CONFIG=
HEPMC2_BIN=/user/local/
HEPMC2_INCLUDE=-I/user/local
HEPMC2_LIB=-L/user/local -Wl,-rpath,/user/local -lHepMC

HEPMC3_USE=false
HEPMC3_CONFIG=
HEPMC3_BIN=
HEPMC3_INCLUDE=
HEPMC3_LIB=

LHAPDF5_USE=false
LHAPDF5_CONFIG=
LHAPDF5_BIN=
LHAPDF5_INCLUDE=
LHAPDF5_LIB=

LHAPDF6_USE=true
LHAPDF6_CONFIG=lhapdf-config
LHAPDF6_BIN=/user/local/
LHAPDF6_INCLUDE=-I/user/local
LHAPDF6_LIB=-L/user/local -Wl,-rpath,/user/local -lLHAPDF

POWHEG_USE=false
POWHEG_CONFIG=lhapdf-config
POWHEG_BIN=
POWHEG_INCLUDE=
POWHEG_LIB=

RIVET_USE=false
RIVET_CONFIG=lhapdf-config
RIVET_BIN=
RIVET_INCLUDE=
RIVET_LIB=

ROOT_USE=true
ROOT_CONFIG=root-config
ROOT_BIN=/home/suman/HEP_Package/ROOT616build/bin/
ROOT_INCLUDE=-I/home/suman/HEP_Package/ROOT616build/include
ROOT_LIB=-L/home/suman/HEP_Package/ROOT616build/lib -Wl,-rpath,/home/suman/HEP_Package/ROOT616build/lib -lCore

YODA_USE=false
YODA_CONFIG=root-config
YODA_BIN=
YODA_INCLUDE=
YODA_LIB=

GZIP_USE=false
GZIP_CONFIG=root-config
GZIP_BIN=
GZIP_INCLUDE=
GZIP_LIB=

PYTHON_USE=false
PYTHON_CONFIG=root-config
PYTHON_BIN=
PYTHON_INCLUDE=
PYTHON_LIB=

MG5MES_USE=false
MG5MES_CONFIG=root-config
MG5MES_BIN=
MG5MES_INCLUDE=
MG5MES_LIB=

OPENMP_USE=false
OPENMP_CONFIG=root-config
OPENMP_BIN=
OPENMP_INCLUDE=
OPENMP_LIB=




#ROOT configuration.
#ROOT_USE=true
#ROOT_USE=false
ROOT_BIN=
ROOT_INCLUDE=./
ROOT_LIB=./

# FASTJET3 configuration.
#FASTJET3_USE=true
FASTJET3_BIN=
FASTJET3_INCLUDE=./
FASTJET3_LIB=./


# GZIP configuration.
#GZIP_USE=false
GZIP_BIN=
GZIP_INCLUDE=./
GZIP_LIB=./

CLHEPSYS=/home/suman/HEP_Package/CLHEPBUILD

CXX=g++ EventShape_vector.cc
CXX_COMMON=-O2  -pedantic -W -Wall -Wshadow -fPIC -lCLHEP-Vector-2.4.1.0
CXX_SHARED=-shared
CXX_SONAME=-Wl,-soname,
LIB_SUFFIX=.so

# Handle GZIP support.
ifeq ($(GZIP_USE),true)
  CXX_COMMON+= -DGZIPSUPPORT -I$(GZIP_INCLUDE) -I${CLHEPSYS} 
  CXX_COMMON+= -L$(GZIP_LIB) -Wl,-rpath,$(GZIP_LIB) -lz -L${CLHEPSYS}
	@echo "I am Inside"
endif

# Check distribution (use local version first, then installed version).
#ifneq ("$(wildcard ../lib/libpythia8.*)","")
#  PREFIX_LIB=../lib
#  PREFIX_INCLUDE=../include
#endif

CXX_COMMON:=-I$(PREFIX_INCLUDE) $(CXX_COMMON)
CXX_COMMON+= -L$(PREFIX_LIB) -Wl,-rpath,$(PREFIX_LIB) -lpythia8 -ldl


HardQCD: HardQCD.cc $(PREFIX_LIB)/libpythia8.a
	$(CXX) $< -o $@ -w -I$(FASTJET3_INCLUDE) $(CXX_COMMON)\
         -L$(FASTJET3_LIB) -Wl,-rpath,$(FASTJET3_LIB) -lfastjet -I$(ROOT_INCLUDE) -I$/usr/local/lib\
         `$(ROOTBIN)root-config --cflags`\
         -Wl,-rpath,$(ROOT_LIB) `$(ROOT_BIN)root-config --glibs`





