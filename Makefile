
FDDIR := $(FDMODULE)
LALIB := lapack
# LALIB := openblas

CXXFLAGS := -O2 -g -Wall -fPIC -Wno-maybe-uninitialized
CXXFLAGS += $(shell root-config --cflags)

LDFLAGS := -O2 -g
LDFLAGS += $(shell root-config --libs) -lMinuit -lGeom

TARG := test
OBJS := Fumili mconvd mtrx_inv_$(LALIB)

run: $(TARG)
	@./$<

test_inv: $(addsuffix .o,test_inv mtrx_inv_$(LALIB) test_mat)
	@echo 'Linking executable $@'
	@$(CXX) $^ $(LDFLAGS) -l$(LALIB) -o $@
	@./$@

%: %.cc
%: %.o
%: %.o $(addsuffix .o,$(OBJS))
	@echo 'Linking executable $@'
	@$(CXX) $^ $(LDFLAGS) -l$(LALIB) -o $@

%.o: %.cc
%.o: %.cc %.d
	@echo 'Compiling $@'
	@$(CXX) $< -c $(CXXFLAGS) $(CPPFLAGS) -o $@

%.d: %.cc
	@echo Making dependency for file $< ...
	@$(CXX) $< -MT $(subst .d,.o,$@) -MM $(CXXFLAGS) -MF $@

clean:
	@echo 'Cleaning'
	@$(RM) $(TARG) test_inv *.o *.d

.PRECIOUS: %.d %.o

.PHONY: run clean

include $(wildcard *.d)
