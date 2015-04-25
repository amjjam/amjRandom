includedir=$(prefix)/include

build:

install: $(includedir)/Random.H

$(includedir)/Random.H: Random.H
	cp Random.H $(includedir)

uninstall:
	- rm -f $(includedir)/Random.H

clean:
