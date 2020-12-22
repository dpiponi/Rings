test: test.cpp
	    c++ -std=c++17 -O3 -I/usr/local/include -L/usr/local/lib -lraylib -F/Library/Developer/CommandLineTools/SDKs/MacOSX10.15.sdk/System/Library/Frameworks -framework OpenGL -framework Cocoa -framework IOKit -framework CoreFoundation -framework CoreVideo -o test test.cpp
