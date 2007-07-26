INCLUDES =  -I$(top_builddir) -I$(top_srcdir)
AM_YFLAGS = -d
AM_LFLAGS = -CF
BUILT_SOURCES = seqparse.h matparse.h

bin_PROGRAMS = ngila

ngila_LDADD = $(LIBOBJS) @LEXLIB@

CLEANFILES = matparse.cc matparse.h matscan.cc seqparse.cc seqparse.h seqscan.cc
noinst_HEADERS = ngila.cmds ngila.h seqdb.h seqparser.h
ngila_SOURCES = align.cpp ngila.cpp seqdb.cpp seqparser.cpp
