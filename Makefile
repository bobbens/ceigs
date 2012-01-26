
PATH_INCLUDE   := /usr/include
PATH_INSTALL   := /usr/lib

LIBNAME	:= libceigs
VERSION	:= 1.1

USE_UMFPACK := true

CFLAGS   := -O3 -fPIC
CFLAGS   += -W -Wall -Wextra -Wunused -Wshadow -Wpointer-arith -Wmissing-prototypes -Winline -Wcast-align -Wmissing-declarations -Wredundant-decls -Wno-long-long -Wcast-align
CFLAGS   += -I/usr/include/suitesparse
LDFLAGS	:= -larpack -lcxsparse
ifeq ($(USE_UMFPACK),true)
CFLAGS	+= -DUSE_UMFPACK
LDFLAGS  += -lumfpack
endif

OBJS		:= ceigs.o ceigs_cs.o cs_fact.o ceigs_lu.o ceigs_cholesky.o ceigs_qr.o ceigs_umfpack.o

.PHONY: all $(LIBNAME) install clean uninstall docs help

all: $(LIBNAME) test

help:
	@echo "Valid targets are:"
	@echo "          all - Makes the library and test application"
	@echo "     libceigs - Makes the ceigs library"
	@echo "         test - Compiles the test application"
	@echo "      install - Installs the library"
	@echo "    uninstall - Uninstalls the library"
	@echo "         docs - Makes the documentation"
	@echo "        clean - Cleans up the build system"
	@echo "         help - Displays this message"

$(LIBNAME): $(LIBNAME).a $(LIBNAME).so

$(LIBNAME).a: $(OBJS)
	$(AR) rcs $(LIBNAME).a $(OBJS)

$(LIBNAME).so: $(OBJS)
	$(CC) -shared -Wl,-soname,$(LIBNAME).so -o $(LIBNAME).so.$(VERSION) $(OBJS)
	ln -sf $(LIBNAME).so.$(VERSION) $(LIBNAME).so

test: $(OBJS) test.o
	$(CC) $(CFLAGS) $(LDFLAGS) $(OBJS) test.o -o test
	./test

install: all
	install -m 644 $(LIBNAME).a $(PATH_INSTALL)
	test -d $(PATH_INCLUDE) || mkdir $(PATH_INCLUDE)
	cp ceigs.h $(PATH_INCLUDE)/ceigs.h
	cp $(LIBNAME).so.$(VERSION) $(PATH_INSTALL)
	(cd $(PATH_INSTALL); ln -sf $(PATH_INSTALL)/$(LIBNAME).so.$(VERSION) $(LIBNAME).so)
	ldconfig

uninstall:
	$(RM) $(PATH_INCLUDE)/ceigs.h
	$(RM) $(PATH_INSTALL)/$(LIBNAME).so.$(VERSION)
	$(RM) $(PATH_INSTALL)/$(LIBNAME).so
	$(RM) $(PATH_INSTALL)/$(LIBNAME).a

docs:
	doxygen
	$(MAKE) -C docs/latex
	cp docs/latex/refman.pdf $(LIBNAME).pdf

clean:
	$(RM) $(OBJS) $(LIBNAME).so.$(VERSION) $(LIBNAME).so $(LIBNAME).a test test.o

