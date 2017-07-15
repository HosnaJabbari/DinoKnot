ACLOCAL_AMFLAGS  = -I m4 --install

bin_PROGRAMS= HFold_interacting_multimodel

AM_CFLAGS  = -g -I./simfold/include -I./simfold/src/common -I./simfold/src/simfold -I ./ -I/usr/include/boost -w -std=c++11

LDFLAGS_BOOST = \
    $(BOOST_SYSTEM_LDFLAGS) \
    $(BOOST_PROGRAM_OPTIONS_LDFLAGS)

LIBS_BOOST = \
    $(BOOST_SYSTEM_LIBS) \
    $(BOOST_PROGRAM_OPTIONS_LIBS) \
    $(BOOST_FILESYSTEM_LIBS)


AM_LDFLAGS  = -g -L./ -lsimfold -lm -lpthread $(LDFLAGS_BOOST) -lboost_regex

# Sources to compile HFold_interacting_multimodel from
HFold_interacting_multimodel_SOURCES = \
	W_final.h W_final.cpp \
	pseudo_loop.h pseudo_loop.cpp \
	V_final.h V_final.cpp \
	VM_final.h VM_final.cpp \
	h_common.h h_common.cpp \
	hfold_validation.h hfold_validation.cpp \
	hfold_interacting.h HFold_interacting.cpp
