CXX = g++
CC = gcc -fms-extensions 
CFLAGS = -Wall -g -O0
INCLUDE = -Iinclude

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	LIBS =  -lm
	OPENCL = -I/System/Library/Frameworks/OpenCL.framework/Versions/A/Headers -L/System/Library/Frameworks/OpenCL.framework/Versions/A/Libraries -framework OpenCL -D_POSIX_C_SOURCE=199309
else
	LIBS =  -lm -lrt
	OPENCL = -I$(OPENCL_ROOT)/include -L$(OPENCL_ROOT)/lib/x86_64 -lOpenCL -D_POSIX_C_SOURCE=199309
endif

all: sampo

lib_icl: src/lib_icl_ext.c src/lib_icl.c
	$(CC) $(CFLAGS) src/lib_icl.c  $(INCLUDE) $(OPENCL) -std=c99 -c -o bin/lib_icl.o 
	$(CC) $(CFLAGS) src/lib_icl_ext.c $(INCLUDE) $(OPENCL) -std=c99 -c -o bin/lib_icl_ext.o 

sampo: src/abms.c lib_icl
	$(CC) $(CFLAGS) src/abms.c src/boltScan.c bin/lib_icl.o bin/lib_icl_ext.o -std=c99 $(INCLUDE) $(LIBS) $(OPENCL) -o bin/sampo

clean:
	rm bin/abms bin/*.o
