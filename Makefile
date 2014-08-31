include make.inc

STINGER_SRC= src/stinger.c src/stinger-utils.c src/timer.c src/xmalloc.c src/stinger-physmap.c src/stinger-iterator.c src/stinger-deprecated.c
CFLAGS+=-Iinclude

ifeq ($(shell uname -s),Linux)
  PYTHON_LIB= -I/usr/include/python2.7
  PYTHON_LIB_FLAG=
  JAVA_LIB= -I /usr/lib/jvm/default-java/include/
  JAVA_LIBNAME=libstinger.so
else 
ifeq ($(shell uname -s),Darwin)
  PYTHON_LIB= -I/System/Library/Frameworks/Python.framework/Headers/
  PYTHON_LIB_FLAG=-lpython
  JAVA_LIB= -I/System/Library/Frameworks/JavaVM.framework/Headers
  JAVA_LIBNAME=libstinger.jnilib
else
  $(warning OS Not supported for java library)
  $(warning OS Not supported for python library)
endif
endif

EXAMPLES=connected-components-static breadth-first-search streaming-clustering-coefficients

.PHONY:	all
ifdef FOR_XMT
all: main gen-streams examples 
else
all: main gen-streams examples rmatter lib
endif

examples: $(EXAMPLES)

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
else
  STINGER_SRC+=obj/x86-full-empty.o
endif

main:	main.c $(STINGER_SRC) $(BLECHIO)
	$(CC) $(MAINPLFLAG) $(CPPFLAGS) $(CFLAGS)  -o $@ $^ \
		$(LDFLAGS) $(LDLIBS)

lib: 	$(subst src,obj,$(subst .c,.o,$(STINGER_SRC))) obj/x86-full-empty.o
	ar rcs libstinger.a obj/*.o

rmatter: src/rmatter.c src/random.c src/timer.c
	$(CC) $(MAINPLFLAG) $(CPPFLAGS) $(CFLAGS) -o $@ $^ \
		$(LDFLAGS) $(LDLIBS)

obj/x86-full-empty.o:	src/x86-full-empty.c
	$(CC) $(MAINPLFLAG) $(CPPFLAGS) $(subst O1,O0,$(subst O2,O0,$(subst O3,O0,$(CFLAGS)))) -c -o obj/x86-full-empty.o $^

static_example:	main.c r-utils.c src/timer.c src/xmalloc.c src/stinger-iterator.c $(BLECHIO) obj/x86-full-empty.o
	$(CC)  $(CFLAGS) -o static_main main.c libstinger.a 


gen-streams:	genstreams/gen-streams.c genstreams/prng.c genstreams/rmat.c src/timer.c src/xmalloc.c $(BLECHIOGEN) $(GETOPT_SRC)
	$(CC) $(GENSTREAMSPLFLAG) $(CPPFLAGS) $(CFLAGS) -o $@ $^ \
		$(LDFLAGS) $(LDLIBS)

xmt-luc-blech.o:	compat/xmt-luc-blech.cc
	$(CXX) $(MAINPLFLAG) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<
xmt-luc-blech-gen.o:	compat/xmt-luc-blech.cc
	$(CXX) $(GENSTREAMPLFLAG) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

connected-components-static:	examples/connected_components_static/connected-comp.c $(STINGER_SRC) $(BLECHIO)
	$(CC) $(MAINPLFLAG) $(CPPFLAGS) $(CFLAGS) -o $@ $^ \
		$(LDFLAGS) $(LDLIBS)

breadth-first-search:	examples/breadth_first_search/bfs.c $(STINGER_SRC) $(BLECHIO)
	$(CC) $(MAINPLFLAG) $(CPPFLAGS) $(CFLAGS) -o $@ $^ \
		$(LDFLAGS) $(LDLIBS)

streaming-clustering-coefficients:	examples/clustering_coefficients_streaming/coeff-main.c examples/clustering_coefficients_streaming/bloom.c examples/clustering_coefficients_streaming/clustering-coeff.c examples/clustering_coefficients_streaming/simple-update.c $(STINGER_SRC) $(BLECHIO)
	$(CC) $(MAINPLFLAG) $(CPPFLAGS) $(CFLAGS) -o $@ $^ \
		$(LDFLAGS) $(LDLIBS)

test_physmap: src/stinger-physmap.c
	$(CC) $(MAINPLFLAG) -D PHYSMAP_TEST -Iinclude -o $@ $^

python: $(subst src,obj,$(subst .c,.o,$(STINGER_SRC))) obj/x86-full-empty.o
	doxygen
	touch python/doxy/doxy.i; \
	rm python/doxy/doxy.i; \
	for f in doc/xml/*.xml; do \
	  python python/doxy2swig.py $$f `echo $$f | sed -e 's/doc\/xml/python\/doxy/' | sed -e 's/.xml/.i/'`; \
	  echo '%include "'`echo $$f | sed -e 's/doc\/xml/doxy/' | sed -e 's/.xml/.i/'`'"' >> python/doxy/doxy.i; \
	done
	touch make.inc; rm make.inc; ln -s make.inc-gcc-openmp make.inc
	swig -python -outdir python -Iinclude -o python/python_stinger_wrapper.c python/stinger.i
	$(CC) -c $(CFLAGS) -fpic -o obj/python_stinger_wrapper.o python/python_stinger_wrapper.c $(PYTHON_LIB) $(PYTHON_LIB_FLAG)
	g++ -shared -fopenmp $^ obj/python_stinger_wrapper.o -o python/_stinger.so $(PYTHON_LIB_FLAG)

java: $(subst src,obj,$(subst .c,.o,$(STINGER_SRC))) obj/x86-full-empty.o
	doxygen
	touch java/doxy/doxy.i; \
	rm java/doxy/doxy.i; \
	for f in doc/xml/*.xml; do \
	  python java/doxy2swig.py $$f `echo $$f | sed -e 's/doc\/xml/java\/doxy/' | sed -e 's/.xml/.i/'`; \
	  echo '%include "'`echo $$f | sed -e 's/doc\/xml/doxy/' | sed -e 's/.xml/.i/'`'"' >> java/doxy/doxy.i; \
	  sed -ibak 's/^\(%[^"]*\)\"/\1 \"\/** \n* /g' `echo $$f | sed -e 's/doc\/xml/java\/doxy/' | sed -e 's/.xml/.i/'`; \
	  sed -ibak 's/\";/\n*\/ public \"/' `echo $$f | sed -e 's/doc\/xml/java\/doxy/' | sed -e 's/.xml/.i/'`; \
	  sed -ibak 's/^\([^/%*]\)/* \1/' `echo $$f | sed -e 's/doc\/xml/java\/doxy/' | sed -e 's/.xml/.i/'`; \
	done
	swig -java -package stinger -outdir java/stinger -Iinclude -o java/java_stinger_wrapper.c java/stinger.i
	$(CC) -c $(CFLAGS) -fpic -o obj/java_stinger_wrapper.o java/java_stinger_wrapper.c $(JAVA_LIB)
	g++ -shared -fopenmp $^ obj/java_stinger_wrapper.o -o java/$(JAVA_LIBNAME)

obj/%.o: src/%.c
	$(CC) $(MAINPLFLAG) $(CPPFLAGS) $(CFLAGS) -fPIC -c $^ -o $@ \
		$(LDFLAGS) $(LDLIBS)

.PHONY: release
release:
	svn co http://svn.cc.gatech.edu/graphs/stinger/trunk trunk
	echo stinger-r`svn info trunk | grep "Revision: " | sed 's/.*Revision: //'` > RELEASE_DIR
	mv trunk `cat RELEASE_DIR`
	tar -czf `cat RELEASE_DIR`.tar.gz `cat RELEASE_DIR`
	rm -rf `cat RELEASE_DIR` RELEASE_DIR

.PHONY:	clean
clean:
	rm -f main gen-streams rmatter $(EXAMPLES) $(BLECHIO) $(BLECHIOGEN) \
		$(MAINPL) $(GENSTREAMSPL) libstinger.a obj/*.o

