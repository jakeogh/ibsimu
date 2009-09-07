SUBDIRS = src tests doc
pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = ibsimu.pc

ACLOCAL_AMFLAGS = -I m4

.PHONY: doc
doc:
	$(MAKE) -C doc doc

CLEANFILES = *~
