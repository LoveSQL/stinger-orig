include ../make.inc

STINGER_SRC= ../src/stinger.c ../src/stinger-utils.c ../src/timer.c ../src/xmalloc.c ../src/stinger-deprecated.c
LIBSTINGER=libstinger.a
LIBCHECK_DIR=check-0.9.8
LIBCHECK=$(LIBCHECK_DIR)/lib/libcheck.a
TEST_BIN=stinger-basic-tests
CFLAGS+=-Iinclude -I../include -I$(LIBCHECK_DIR)/include


.PHONY:	all
all: test

ifdef USE_GETOPT_NETBSD
  GETOPT_SRC=genstreams/getopt-netbsd/getopt_long.c
  CPPFLAGS+=-Igenstreams/getopt-netbsd
endif

ifdef FOR_XMT
  BLECHIO=xmt-luc-blech.o
  BLECHIOGEN=xmt-luc-blech-gen.o
  MAINPL=main.pl
  MAINPLFLAG=-pl $(MAINPL)
  GENSTREAMSPL=gen-streams.pl
  GENSTREAMSPLFLAG=-pl $(GENSTREAMSPL)
endif

libcheck.a:
	cd $(LIBCHECK_DIR); ./configure --prefix=`pwd`; make install

libstinger.a: 	$(STINGER_SRC) ../obj/x86-full-empty.o
	$(CC) -c $(CFLAGS) $(STINGER_SRC) -idirafter "../include" 
	ar -cvq libstinger.a stinger.o stinger-utils.o timer.o xmalloc.o ../obj/x86-full-empty.o stinger-deprecated.o
	rm *.o

test:	libcheck.a libstinger.a basic-tests.cc
	$(CXX) $(CFLAGS) -Drestrict=__restrict__ basic-tests.cc -o $(TEST_BIN) $(LIBSTINGER) $(LIBCHECK) $(LDLIBS)

../obj/x86-full-empty.o:	../src/x86-full-empty.c
	$(CC) -I../include -c -o ../obj/x86-full-empty.o $^ $(LDFLAGS) $(LDLIBS)

.PHONY:	clean distclean
clean:
	rm -f $(TEST_BIN) *.o $(LIBSTINGER) $(MAINPL)
	-if test -f $(LIBCHECK_DIR)/Makefile ; then make -C $(LIBCHECK_DIR) clean; fi

distclean:
	rm -f $(TEST_BIN) *.o $(LIBSTINGER) $(MAINPL)
	-if test -f $(LIBCHECK_DIR)/Makefile ; then make -C $(LIBCHECK_DIR) distclean; fi
