
OFILES=mpfit.o pplm_mpfit.o
LIBFILE=libmpfit.a
CFLAGS+=-fopenmp

MPFIT_FILES   = mpfit.c mpfit.h testmpfit.c
PPLM_FILES    = pplm_mpfit.h pplm_mpfit.c
ADD2TAR_FILES = Makefile testpplm_mpfit.c notes.txt cmpfit-1.4.zip
TAR_FILE      = pplm_mpfit.tgz

all:  cmpfit-1.4 $(LIBFILE) testmpfit testpplm_mpfit

clean:
	rm -f $(OFILES) testmpfit testpplm_mpfit $(LIBFILE)

cleanall: clean
	rm -f $(TAR_FILE)
	rm -f $(MPFIT_FILES)
	rm -rdf cmpfit-1.4

tarball:
	tar -zcvf $(TAR_FILE)  $(PPLM_FILES) $(ADD2TAR_FILES)

cmpfit-1.4:
	unzip $(@).zip

mpfit.o: mpfit.c mpfit.h
	$(CC) -c -o $@ $< $(CFLAGS)

mpfit.c:
	cp cmpfit-1.4/$(@) .

mpfit.h:
	cp cmpfit-1.4/$(@) .

testmpfit.c:
	cp cmpfit-1.4/$(@) .

pplm_mpfit.o: pplm_mpfit.c pplm_mpfit.h
	$(CC) -c -o $@ $< $(CFLAGS)

$(LIBFILE): $(OFILES)
	$(AR) r $@ $(OFILES)

testmpfit: testmpfit.c libmpfit.a
	$(CC) -o $@ $(CFLAGS) testmpfit.c -L. -lmpfit -lm

testpplm_mpfit: testpplm_mpfit.c libmpfit.a
	$(CC) -o $@ $(CFLAGS) testpplm_mpfit.c -L. -lmpfit -lm


