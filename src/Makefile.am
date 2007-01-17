INCLUDES =  -I$(top_builddir) -I$(top_srcdir)
AM_YFLAGS = -d
AM_LFLAGS = -CF
BUILT_SOURCES = seqparse.h matparse.h

bin_PROGRAMS = ngila ngilag

ngila_SOURCES = ngila.cc align.cc seqparse.yy seqscan.ll matparse.yy matscan.ll
ngila_LDADD = $(LIBOBJS) @LEXLIB@

ngilag_SOURCES = ngila.cc aligng.cc seqparse.yy seqscan.ll matparse.yy matscan.ll
ngilag_LDADD = $(LIBOBJS) @LEXLIB@


CLEANFILES = matparse.cc matparse.h matscan.cc seqparse.cc seqparse.h seqscan.cc
noinst_HEADERS = ngila.h

