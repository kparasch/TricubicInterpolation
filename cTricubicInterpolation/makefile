objdir = object
srcdir = source
incdir = include

python3inc = `python3 -c "from distutils import sysconfig; print(sysconfig.get_python_inc())"`
numpy3inc = `python3 -c "import numpy; print(numpy.get_include())"`
python2inc = `python2 -c "from distutils import sysconfig; print(sysconfig.get_python_inc())"`
numpy2inc = `python2 -c "import numpy; print(numpy.get_include())"`

objects = $(objdir)/coefs.o $(objdir)/derivs.o $(objdir)/evaluate.o   

cc = gcc
cflags = -O3 -std=c99 -Wall -fPIC -I$(incdir) 
ldflags =  -lm #-lpython

all: Tricubic_c.so Tricubic2_c.so
py3: Tricubic_c.so
py2: Tricubic_c.so

Tricubic_c.so: $(objects) $(objdir)/pywrap3.o
	@$(cc) --shared -fPIC $(cflags) $^ -o $@ $(ldflags) 

Tricubic_c.so: $(objects) $(objdir)/pywrap2.o
	@$(cc) --shared -fPIC $(cflags) $^ -o $@ $(ldflags) 

$(objdir)/coefs.o: $(srcdir)/coefs.c $(incdir)/coefs.h 
	@$(cc) -c $(cflags) $< -o $@ 

$(objdir)/derivs.o: $(srcdir)/derivs.c $(incdir)/derivs.h 
	@$(cc) -c $(cflags) $< -o $@ 

$(objdir)/evaluate.o: $(srcdir)/evaluate.c $(incdir)/evaluate.h 
	@$(cc) -c $(cflags) $< -o $@ 

$(objdir)/pywrap3.o: $(srcdir)/pywrap.c $(incdir)/pywrap.h 
	@echo python3 include path: $(python3inc)
	@echo numpy \(for python3\) include path: $(numpy3inc)
	@$(cc) -c $(cflags) -I$(python3inc) -I$(numpy3inc) $< -o $@ 

$(objdir)/pywrap2.o: $(srcdir)/pywrap.c $(incdir)/pywrap.h 
	@echo python2 include path: $(python2inc)
	@echo numpy \(for python2\) include path: $(numpy2inc)
	@$(cc) -c $(cflags) -I$(python2inc) -I$(numpy2inc) $< -o $@ 

clean:
	rm $(objdir)/*
	rm Tricubic_c.so
	rm Tricubic2_c.so

