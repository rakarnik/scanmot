#
# Macros
#

CC = /usr/bin/g++
CC_OPTIONS = -O3 -Wall -Wextra
#CC_OPTIONS = -O0 -ggdb -Wall -Wextra
LNK_OPTIONS =

#
# INCLUDE directories
#
INCLUDE = -I. -Isrc


#
# Build scanmot
#

all: scanmot scantop

scanmot: \
		bin/bgmodel.o\
		bin/motif.o\
		bin/scanmot.o\
		bin/scoredsite.o\
		bin/seqset.o\
		bin/site.o\
		bin/siteheap.o\
		bin/standard.o
	$(CC) $(LNK_OPTIONS) \
		bin/bgmodel.o\
		bin/motif.o\
		bin/scanmot.o\
		bin/scoredsite.o\
		bin/seqset.o\
		bin/site.o\
		bin/siteheap.o\
		bin/standard.o\
		-o bin/scanmot

scantop: \
		bin/bgmodel.o\
		bin/motif.o\
		bin/scantop.o\
		bin/seqset.o\
		bin/site.o\
		bin/standard.o
	$(CC) $(LNK_OPTIONS) \
		bin/bgmodel.o\
		bin/motif.o\
		bin/scantop.o\
		bin/seqset.o\
		bin/site.o\
		bin/standard.o\
		-o bin/scantop

clean: 
	rm -f bin/*.o bin/scanmot bin/scantop

bin/bgmodel.o: src/bgmodel.h src/bgmodel.cpp
	$(CC) $(CC_OPTIONS) src/bgmodel.cpp -c $(INCLUDE) -o bin/bgmodel.o

bin/motif.o: src/motif.h src/motif.cpp
	$(CC) $(CC_OPTIONS) src/motif.cpp -c $(INCLUDE) -o bin/motif.o

bin/scanmot.o: src/scanmot.h src/scanmot.cpp
	$(CC) $(CC_OPTIONS) src/scanmot.cpp -c $(INCLUDE) -o bin/scanmot.o

bin/scantop.o: src/scantop.h src/scantop.cpp
	$(CC) $(CC_OPTIONS) src/scantop.cpp -c $(INCLUDE) -o bin/scantop.o

bin/scoredsite.o: src/scoredsite.cpp src/scoredsite.h
	$(CC) $(CC_OPTIONS) src/scoredsite.cpp -c $(INCLUDE) -o bin/scoredsite.o

bin/seqset.o: src/seqset.h src/seqset.cpp
	$(CC) $(CC_OPTIONS) src/seqset.cpp -c $(INCLUDE) -o bin/seqset.o

bin/site.o: src/site.h src/site.cpp
	$(CC) $(CC_OPTIONS) src/site.cpp -c $(INCLUDE) -o bin/site.o

bin/siteheap.o: src/siteheap.h src/siteheap.cpp
	$(CC) $(CC_OPTIONS) src/siteheap.cpp -c $(INCLUDE) -o bin/siteheap.o

bin/standard.o: src/standard.h src/standard.cpp
	$(CC) $(CC_OPTIONS) src/standard.cpp -c $(INCLUDE) -o bin/standard.o
