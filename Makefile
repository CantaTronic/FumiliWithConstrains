
FDDIR := $(FDMODULE)

CXXFLAGS := -O2 -g -Wall -fPIC -Wno-maybe-uninitialized
CXXFLAGS += $(shell root-config --cflags)

LDFLAGS := -O2 -g
LDFLAGS += $(shell root-config --libs) -lMinuit -lGeom

TARG := test
OBJS := Fumili mconvd

run: $(TARG)
	@./$<

%: %.cc
%: %.o
%: %.o $(addsuffix .o,$(OBJS))
	@echo 'Linking executable $@'
	@$(CXX) $^ $(LDFLAGS) -o $@

%.o: %.cc
%.o: %.cc %.d
	@echo 'Compiling $@'
	@$(CXX) $< -c $(CXXFLAGS) $(CPPFLAGS) -o $@

%.d: %.cc
	@echo Making dependency for file $< ...
	@$(CXX) $< -MT $(subst .d,.o,$@) -MM $(CXXFLAGS) -MF $@

clean:
	@echo 'Cleaning'
	@$(RM) $(TARG) $(addsuffix .o,$(TARG)) $(addsuffix .d,$(TARG)) \
	$(addsuffix .o,$(OBJS)) $(addsuffix .d,$(OBJS))

.PRECIOUS: %.d %.o

.PHONY: run clean

include $(wildcard *.d)
