# what flags you want to pass to the C compiler & linker
AM_CXXFLAGS = -Wall -Wnon-virtual-dtor -Wextra -Wpedantic -g -O2 -std=c++17 -pthread

AM_CCFLAGS = -O2 -g

# this lists the binaries to produce, the (non-PHONY, binary) targets in
# the previous manual Makefile
if WEIGHT128
bin_PROGRAMS = bin/kira128


bin_kira128_SOURCES = \
kira/main.cpp \
kira/connect2kira.cpp \
kira/dataBase.cpp \
kira/dgl.cpp \
kira/ibp.cpp \
kira/integral.cpp \
kira/kira.cpp \
kira/ReadYamlFiles.cpp \
kira/reduction.cpp \
kira/threadpool.cpp \
kira/tools.cpp \
kira/trivial_sym.cpp \
pyred/coeff_helper.cpp \
pyred/coeff_int.cpp \
pyred/coeff_vec.cpp \
pyred/integrals.cpp \
pyred/interface.cpp \
pyred/keyvaluedb.cpp \
pyred/parser.cpp \
pyred/relations.cpp \
sqlite3/sqlite3.c \
gzstream/gzstream.cpp \
kira/black_nox.h \
kira/connect2kira.h \
kira/createSeeds.h \
kira/dataBase.h \
kira/integral.h \
kira/kira.h \
kira/ReadYamlFiles.h \
kira/threadpool.h \
kira/tools.h \
kira/trivial_sym.h \
pyred/coeff_helper.h \
pyred/coeff_int.h \
pyred/coeff_vec.h \
pyred/defs.h \
pyred/gauss.h \
pyred/integrals.h \
pyred/interface.h \
pyred/keyvaluedb.h \
pyred/parser.h \
pyred/relations.h \
sqlite3/sqlite3.h \
gzstream/gzstream.h

if FIREFLY
bin_kira128_SOURCES += kira/reductionFF.cpp \
kira/mpi_work.cpp
endif

bin_kira128_CPPFLAGS = -DGZSTREAM_NAMESPACE=kira_gzstream $(FIREFLY_CFLAGS) $(GINAC_CFLAGS) $(CLN_CFLAGS) $(YAML_CPP_CFLAGS) $(DL_CFLAGS)

if FIREFLY
bin_kira128_CPPFLAGS += -DKIRAFIREFLY
endif

if DEBUG
bin_kira128_CPPFLAGS += -DVALGRIND
endif

if WEIGHT128
bin_kira128_CPPFLAGS += -DWEIGHT128
endif

bin_kira128_LDADD = $(FIREFLY_LIBS) $(GINAC_LIBS) $(CLN_LIBS) $(YAML_CPP_LIBS) $(DL_LIBS) -lstdc++ -L. -lz

bin_kira128_LDFLAGS =
else
bin_PROGRAMS = bin/kira


bin_kira_SOURCES = \
kira/main.cpp \
kira/connect2kira.cpp \
kira/dataBase.cpp \
kira/dgl.cpp \
kira/ibp.cpp \
kira/integral.cpp \
kira/kira.cpp \
kira/ReadYamlFiles.cpp \
kira/reduction.cpp \
kira/threadpool.cpp \
kira/tools.cpp \
kira/trivial_sym.cpp \
pyred/coeff_helper.cpp \
pyred/coeff_int.cpp \
pyred/coeff_vec.cpp \
pyred/integrals.cpp \
pyred/interface.cpp \
pyred/keyvaluedb.cpp \
pyred/parser.cpp \
pyred/relations.cpp \
sqlite3/sqlite3.c \
gzstream/gzstream.cpp \
kira/black_nox.h \
kira/connect2kira.h \
kira/createSeeds.h \
kira/dataBase.h \
kira/integral.h \
kira/kira.h \
kira/ReadYamlFiles.h \
kira/threadpool.h \
kira/tools.h \
kira/trivial_sym.h \
pyred/coeff_helper.h \
pyred/coeff_int.h \
pyred/coeff_vec.h \
pyred/defs.h \
pyred/gauss.h \
pyred/integrals.h \
pyred/interface.h \
pyred/keyvaluedb.h \
pyred/parser.h \
pyred/relations.h \
sqlite3/sqlite3.h \
gzstream/gzstream.h

if FIREFLY
bin_kira_SOURCES += kira/reductionFF.cpp \
kira/mpi_work.cpp
endif

bin_kira_CPPFLAGS = -DGZSTREAM_NAMESPACE=kira_gzstream $(FIREFLY_CFLAGS) $(GINAC_CFLAGS) $(CLN_CFLAGS) $(YAML_CPP_CFLAGS) $(DL_CFLAGS)

if FIREFLY
bin_kira_CPPFLAGS += -DKIRAFIREFLY
endif

if DEBUG
bin_kira_CPPFLAGS += -DVALGRIND
endif

bin_kira_LDADD = $(FIREFLY_LIBS) $(GINAC_LIBS) $(CLN_LIBS) $(YAML_CPP_LIBS) $(DL_LIBS) -lstdc++ -L. -lz

bin_kira_LDFLAGS =
endif
