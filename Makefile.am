SUBDIRS = src tests doc
pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = ibsimu.pc

ACLOCAL_AMFLAGS = -I m4

.PHONY: doc webdoc

doc:
	$(MAKE) -C doc doc

webdoc:
	$(MAKE) -C doc webdoc


EXTRA_DIST = NEWS COMPATIBILITY

CLEANFILES = *~
