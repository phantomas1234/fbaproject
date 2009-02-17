swig -python glpki.i;
gcc -DNDEBUG -g -O3 -I/usr/local/include -I/Users/niko/arbeit/Software/glpk-4.24/include -I/Library/Frameworks/Python.framework/Versions/2.5/include/python2.5 -c glpki_wrap.c -o glpki_wrap.o;
gcc -arch i386 -isysroot /Developer/SDKs/MacOSX10.4u.sdk -g -bundle -undefined dynamic_lookup -L/usr/local/lib -I/usr/local/include glpki_wrap.o -L/usr/local/lib -lglpk -o _glpki.so