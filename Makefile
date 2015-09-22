#
# Macros
#

CC = /usr/bin/g++
CC_OPTIONS = -O3 -Wall -Wextra
LNK_OPTIONS =
BIN_DIR = bin

#
# INCLUDE directories
#
INCLUDE = -I. -Isrc


#
# Build scanmot
#

all: scanmot

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

clean: 
	rm -f bin/*.o bin/scanmot

dir_guard=@mkdir -p $(@D)

$(BIN_DIR)/%.o: src/%.cpp
	$(dir_guard)
	$(CC) $(CC_OPTIONS) -c $(INCLUDE) -o $@ $<


