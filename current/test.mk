DEPTRACKING=-MD -MF $(@:.o=.d)
CXXFLAGS:=-g -std=c++14 -Wall -Wextra -Isrc
BUILDEXE=g++ -o$@ $(CXXFLAGS) $(LDFLAGS) $^

CHECKDIR=@mkdir -p $(dir $@)

all: bin/%.o


bin/%.o: src/%.cpp
	$(CHECKDIR)
	g++ -o$@ -c $(CXXFLAGS) $(DEPTRACKING) $<

bin/examples/%.o: examples/%.cpp
	$(CHECKDIR)
	g++ -o$@ -c $(CXXFLAGS) $(DEPTRACKING) $<

clean:
	#find bin -name '*.d' -delete -o -name '*.o' -delete -o '(' -perm -u=x '!' -type d ')' -delete