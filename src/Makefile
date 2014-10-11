INSTALLDIR=../..

LIBINSTALLDIR=$(INSTALLDIR)/lib/
INCLUDEINSTALLDIR=$(INSTALLDIR)/include/

RANDOM_SRC=random.C
RANDOM_OBJ=$(RANDOM_SRC:.C=.o)

install: librandom.a
	- mkdir $(LIBINSTALLDIR)
	cp librandom.a $(LIBINSTALLDIR)
	- mkdir $(INCLUDEINSTALLDIR)
	cp ../include/random.h $(INCLUDEINSTALLDIR)

librandom.a: $(RANDOM_OBJ)
	ar rc $@ $^
