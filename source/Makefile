# Generic instructions:  make, then make flca, then make warp,
# then as root, make libflca-install, make flcainstall, and make warpinstall.
# To uninstall, as root, make flcalib-uninstall, then make flcauninstall, then
# make warpuninstall .
# To remove executable, objects, and libflca.a from build directory, make clean.
# You can edit this file to change the choice of compiler CC, 
# install directories
# FLCA_BINDIR, FLCA_MANDIR, FLCA_INCLUDEDIR, FLCA_LIBDIR, compiler 
# options COPTS, or location of
# fftw3 library LIBFFTW3, or location of FFTW include files INCLUDEFFTW3 . 
# If no root permission, will need to change these locations. 
#CC = gcc4
#CC = gcc
CC = gcc
#CC = tcc
# linux:
#FLCA_MANDIR = /usr/local/man/man1
# macosx:
#FLCA_MANDIR = /usr/share/man/man1
FLCA_MANDIR = /usr/local/share/man/man1
FLCA_BINDIR = /usr/local/bin
LIBFFTW3 = /usr/local/lib
FLCA_LIBDIR = /usr/local/lib
FLCA_INCLUDEDIR = /usr/local/include
LOPTS = -lflca -lfftw3 -lm
INCLUDEFFTW3 = /usr/local/include
# next version of options for a reasonably optimized executable
#COPTS = -O3 -Wall -fomit-frame-pointer -fPIC
COPTS = -O3 -Wall -fomit-frame-pointer -fPIC
# next version of options for debugging with gdb or valgrind
#COPTS = -g
libflca.a:	flcasubs.c flcasubs.h
	$(CC) $(COPTS) -shared flcasubs.c -I$(INCLUDEFFTW3) -I. -o flca.so
#	ar -rs libflca.a flcasubs.o
all:	libflca.a
clean:
	rm -f flcasubs.so
