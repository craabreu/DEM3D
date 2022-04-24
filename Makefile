# Define DEBUG or FAST mode:
#   Build the "fast" version with: `make` or `make DEBUG=0`
#   Build the "debug" version with: `make DEBUG=1`
DEBUG ?= 0

# Installation prefix:
PREFIX ?= /usr/local

# Compilers and their basic options:
FORT ?= gfortran
BASIC_F_OPTS = -march=native -cpp -fmax-errors=1 -Wunused

# Option FAST (default):
FAST_F_OPTS = -O3

# Option DEBUG:
DEBUG_F_OPTS = -Wall -Wno-maybe-uninitialized

# Checks chosen option:
ifeq ($(DEBUG), 1)
  FOPTS = $(BASIC_F_OPTS) $(DEBUG_F_OPTS)
else
  FOPTS = $(BASIC_F_OPTS) $(FAST_F_OPTS)
endif

SRCDIR  = ./src
OBJDIR  = $(SRCDIR)/obj
BINDIR  = ./bin

src = $(addprefix $(SRCDIR)/, $(addsuffix .f90, $(1)))
obj = $(addprefix $(OBJDIR)/, $(addsuffix .o, $(1)))
SRC  = $(call src, MTypes MDimension3D MUtil MTime BinFiles MCells MConfig MContacts mFluid MForces MDEM DEM)
OBJ  = $(patsubst $(SRCDIR)/%.f90,$(OBJDIR)/%.o,$(SRC))

all: $(BINDIR)/dem3d

.PHONY: clean

clean:
	rm -rf $(OBJDIR) $(BINDIR)

install:
	cp -f $(BINDIR)/* $(PREFIX)/bin

$(BINDIR)/dem3d: $(OBJ)
	mkdir -p $(BINDIR)
	$(FORT) $(FOPTS) -J$(OBJDIR) -o $@ $^

$(OBJDIR)/DEM.o: $(SRCDIR)/DEM.f90 $(call obj, MDEM MUtil MTime mFluid)
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/MDEM.o: $(SRCDIR)/MDEM.f90 $(call obj, MForces MUtil)
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/MForces.o: $(SRCDIR)/MForces.f90 $(call obj, MTypes MDimension3D MCells MContacts MConfig mFluid)
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/mFluid.o: $(SRCDIR)/mFluid.f90 $(call obj, MTypes MDimension3D MConfig MUtil)
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/MContacts.o: $(SRCDIR)/MContacts.f90 $(call obj, MTypes MDimension3D BinFiles)
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/MConfig.o: $(SRCDIR)/MConfig.f90 $(call obj, MDimension3D MUtil)
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/MCells.o: $(SRCDIR)/MCells.f90 $(call obj, MTypes BinFiles)
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/MUtil.o: $(SRCDIR)/MUtil.f90 $(call obj, MTypes BinFiles)
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/MDimension3D.o: $(SRCDIR)/MDimension3D.f90 $(call obj, MTypes BinFiles)
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/MTime.o: $(SRCDIR)/MTime.f90
	mkdir -p $(OBJDIR)
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/BinFiles.o: $(SRCDIR)/BinFiles.f90
	mkdir -p $(OBJDIR)
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/MTypes.o: $(SRCDIR)/MTypes.f90
	mkdir -p $(OBJDIR)
	$(FORT) $(FOPTS) -J$(OBJDIR) -c -o $@ $<
