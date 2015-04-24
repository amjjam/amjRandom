ifndef prefix
	prefix=~
endif
export prefix

build:
	$(MAKE) -C src build
	$(MAKE) -C include build

install:
	$(MAKE) -C src
	$(MAKE) -C include

clean: 
	$(MAKE) -C src clean
	$(MAKE) -C include clean

uninstall:
	$(MAKE) -C src uninstall
	$(MAKE) -C include uninstall
