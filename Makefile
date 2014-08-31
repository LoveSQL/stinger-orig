include make.inc

STINGER_SRC= src/stinger.c src/stinger-utils.c src/timer.c src/xmalloc.c src/stinger-physmap.c src/stinger-iterator.c src/stinger-deprecated.c
PYTHON_LIB= -I/usr/include/python2.7
JAVA_LIB= -I/usr/lib/jvm/java-7-openjdk-amd64/include/
CFLAGS+=-Iinclude

EXAMPLES=connected-components-static breadth-first-search streaming-clustering-coefficients

.PHONY:	all
all: main gen-streams examples

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

lib: 	$(STINGER_SRC)
	$(CC) -c $(CFLAGS) $(STINGER_SRC) -idirafter "include" 
	ar -cvq libstinger.a stinger.o stinger-utils.o timer.o xmalloc.o
	rm *.o

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
	swig -python -outdir python -Iinclude -o python/python_stinger_wrapper.c python/stinger.i
	$(CC) -c $(CFLAGS) -fpic -o obj/python_stinger_wrapper.o python/python_stinger_wrapper.c $(PYTHON_LIB)
	g++ -shared -fopenmp -lrt $^ obj/python_stinger_wrapper.o -o python/_stinger.so

java: $(subst src,obj,$(subst .c,.o,$(STINGER_SRC))) obj/x86-full-empty.o
	swig -java -outdir java -Iinclude -o java/java_stinger_wrapper.c java/stinger.i
	$(CC) -c $(CFLAGS) -fpic -o obj/java_stinger_wrapper.o java/java_stinger_wrapper.c $(JAVA_LIB)
	g++ -shared -fopenmp -lrt $^ obj/java_stinger_wrapper.o -o java/libstinger.so

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
	rm -f main gen-streams $(EXAMPLES) $(BLECHIO) $(BLECHIOGEN) \
		$(MAINPL) $(GENSTREAMSPL) libstinger.a obj/*.o
