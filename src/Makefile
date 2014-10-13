INSTALLDIR=~/lib/

RANDOM_SRC=random.C
RANDOM_OBJ=$(RANDOM_SRC:.C=.o)

install: librandom.a
	cp librandom.a $(INSTALLDIR)

uninstall:
	rm $(INSTALLDIR)/librandom.a

clean:
	rm $(RANDOM_OBJ)

librandom.a: $(RANDOM_OBJ)
	ar rc $@ $^
