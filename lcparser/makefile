liblcparser.so: lcparserUX.cxx lcparser.cpp
	g++ -shared -fPIC -std=c++17 -o liblcparser.so lcparser.cpp lcparserUX.cxx `root-config --cflags --libs`
	- mkdir lib 2> /dev/null
	- mv liblcparser.so *.pcm lib/ 2> /dev/null

lcparserUX.cxx: lcparser.h linkdef.h
	rootcling -f lcparserUX.cxx -c lcparser.h linkdef.h


clean:
	-rm -f lcparserUX.cxx lcparser.so *.pcm lib/*.pcm lib/liblcparser.so 2> /dev/null