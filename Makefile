
# LALIB := LAPACK
# LALIB := OPENBLAS
LALIB := ROOT

CXXFLAGS := -O2 -g -Wall -fPIC -Wno-maybe-uninitialized
ifeq ($(LALIB),ROOT)
  CXXFLAGS += $(shell root-config --cflags)
endif

LDFLAGS := -O2 -g
ifeq ($(LALIB),LAPACK)
  LDFLAGS += -llapack
endif
ifeq ($(LALIB),OPENBLAS)
  LDFLAGS += -lopenblas
endif
ifeq ($(LALIB),ROOT)
  LDFLAGS += $(shell root-config --libs) -lMinuit -lFumili
endif

TARG := test
OBJS := Fumili mconvd mtrx_inv

run: $(TARG)
	@./$<

$(TARG): $(TARG).o $(addsuffix .o,$(OBJS))
	@echo 'Linking executable $@'
	@$(CXX) $^ $(LDFLAGS) -o $@

TestFitter: TestFitter.o AbstractFitter.o MyPDF.o
	@echo 'Linking executable $@'
	@$(CXX) $^ $(LDFLAGS) -o $@

test_inv: $(addsuffix .o,test_inv mtrx_inv test_mat)
	@echo 'Linking executable $@'
	@$(CXX) $^ $(LDFLAGS) -o $@
	@./$@

unif: $(addsuffix .o, test_unif PDFGen PDFGenSquare MyPDF)
	@echo 'Linking executable $@'
	@$(CXX) $^ $(LDFLAGS) -o $@
# 	@./$@

%: %.cc
%: %.o
%: %.o
	@echo 'Linking executable $@'
	@$(CXX) $^ $(LDFLAGS) -o $@

%.o: %.cc
%.o: %.cc %.d
	@echo 'Compiling $@'
	@$(CXX) $< -c $(CXXFLAGS) -DMTRX_LIB_$(LALIB) -o $@

test_mat.o: test_mat.cc test_mat.d
	@echo 'Compiling $@'
	@$(CXX) $< -c $(CXXFLAGS) -Wno-unused -DMTRX_LIB_$(LALIB) -o $@

%.d: %.cc
	@echo Making dependency for file $< ...
	@$(CXX) $< -MT $(subst .d,.o,$@) -MM $(CXXFLAGS) -DMTRX_LIB_$(LALIB) -MF $@

clean:
	@echo 'Cleaning'
	@$(RM) $(TARG) test_inv *.o *.d

.PRECIOUS: %.d %.o

.PHONY: run clean

include $(wildcard *.d)
