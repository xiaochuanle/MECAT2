PWD := $(shell pwd)
OS_TYPE		:= $(shell echo `uname`)
MACHINE_TYPE	:= $(shell echo `uname -m`)
ifeq (${MACHINE_TYPE}, x86_64)
	MACHINE_TYPE = amd64
endif
BUILD_DIR	:= ${PWD}/${OS_TYPE}-${MACHINE_TYPE}/bin

.PHONY: all clean
all: mecat

mecat:
	cd src && make

clean:
	cd src && make clean
