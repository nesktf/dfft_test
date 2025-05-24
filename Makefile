CXX := g++
CXXFLAGS := -std=c++20
LDLIBS := -lfmt

CXXSRCS := dfft.cpp
CXXOBJS := $(subst .cpp,.o,$(CXXSRCS))
CXXEXE := fft

LUAMAIN := dfft.lua
LUAINT := luajit

cpp: cpp_build
	./$(CXXEXE)

cpp_build: $(CXXOBJS)
	$(CXX) -o $(CXXEXE) $(CXXOBJS) $(LDLIBS)

%.o: %.c
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm $(CXXOBJS)

distclean: clean
	rm $(CXXEXE)

lua: $(LUAMAIN)
	$(LUAINT) $(LUAMAIN)
