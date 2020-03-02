PWD := $(shell pwd)
OS_TYPE		:= $(shell echo `uname`)
MACHINE_TYPE	:= $(shell echo `uname -m`)
ifeq (${MACHINE_TYPE}, x86_64)
	MACHINE_TYPE = amd64
endif
BUILD_DIR	:= ${PWD}/${OS_TYPE}-${MACHINE_TYPE}/bin

BUILD_TOP_DIR := ${PWD}/${OS_TYPE}-${MACHINE_TYPE}

HDF5_VERSION ?= 1.10.4
HDF5_INSTALL = ${BUILD_TOP_DIR}/hdf5
HDF5_SOURCE = ${BUILD_TOP_DIR}/hdf5-${HDF5_VERSION}
HDF5_INCLUDE = ${HDF5_INSTALL}/include
HDF5_LIB = ${HDF5_INSTALL}/lib/libhdf5.a
HDF5_TAR_GZ = ${BUILD_TOP_DIR}/hdf5-$(HDF5_VERSION).tar.gz

DEXTRACTOR_BIN_NAME = dexqv dexta dextract undexqv undexta
DEXTRACTOR_BIN = $(patsubst %, ${BUILD_DIR}/%, $(DEXTRACTOR_BIN_NAME))

.PHONY: all clean mecat dextractor

all: mecat dextractor

mecat:
	cd src && make

clean:
	cd src && make clean
	cd DEXTRACTOR && make -f ../dextract_makefile clean
	rm ${HDF5_INSTALL} -rf
	rm ${HDF5_SOURCE} -rf
	rm ${DEXTRACTOR_BIN} -f

dextractor: ${BUILD_TOP_DIR}/bin/dexqv ${BUILD_TOP_DIR}/bin/dexta ${BUILD_TOP_DIR}/bin/dextract ${BUILD_TOP_DIR}/bin/undexqv ${BUILD_TOP_DIR}/bin/undexta


${DEXTRACTOR_BIN}: ${HDF5_LIB}
	mkdir -p ${BUILD_DIR}		
	cd DEXTRACTOR && make -f ../dextract_makefile HDF5_INCLUDE=${HDF5_INCLUDE} HDF5_LIB=${HDF5_LIB} 
	cd DEXTRACTOR && cp ${DEXTRACTOR_BIN_NAME} ${BUILD_DIR}

${HDF5_LIB}: ${HDF5_TAR_GZ}
	tar -xzf ${BUILD_TOP_DIR}/hdf5-$(HDF5_VERSION).tar.gz -C ${BUILD_TOP_DIR} || exit 255
	cd ${BUILD_TOP_DIR}/hdf5-$(HDF5_VERSION) && \
		./configure --enable-threadsafe --disable-hl --libdir=`pwd`/../hdf5/lib --includedir=`pwd`/../hdf5/include --prefix=`pwd`/../hdf5 && \
		make -j ${MAKEFLAGS} && make install

${HDF5_TAR_GZ}:
	mkdir -p ${BUILD_DIR}		
	version_major_minor=`echo "$(HDF5_VERSION)" | sed -E 's/\.[0-9]+$$//'`; \
	cp third_party/hdf5-1.10.4.tar.gz ${BUILD_TOP_DIR}/
	

