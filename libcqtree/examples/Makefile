include ../Make.helper
CXX_FLAGS=$(MY_CXX_FLAGS) $(MY_CXX_OPT_FLAGS) \
-I$(INC_DIR) \
-I$(INC_DIR)/../../cdswrapper/include/ \
-I$(INC_DIR)/../../../libcds/include/ \
-I$(INC_DIR)/../../../sdsl-lite/include/

CCLIB=$(LIB_DIR)/libcqtree.a \
$(LIB_DIR)/ext_cdswrapper/ext_sdsl-lite/lib/libsdsl.a \
$(LIB_DIR)/ext_cdswrapper/ext_libcds/libcds.a \
$(LIB_DIR)/ext_cdswrapper/libcdswrapper.a 

SOURCES=$(wildcard *.cc)
EXECS=$(SOURCES:.cc=.x)

all: $(EXECS)

build-test: $(EXECS)
	        
%.x:%.cc
	$(MY_CXX) $(CXX_FLAGS) -o $@ $< $(CCLIB) 

clean:
	rm -f $(EXECS)
	rm -rf *.dSYM
