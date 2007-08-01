AM_CPPFLAGS= -I$(top_srcdir)/include $(all_includes)

bin_PROGRAMS = ngila

ngila_CPPFLAGS = $(BOOST_CPPFLAGS) $(GSL_CFLAGS)
ngila_LDFLAGS = $(all_libraries) $(BOOST_LDFLAGS) $(ALL_STATIC_LDFLAGS)
ngila_LDADD = ${LIBOJS} $(BOOST_PROGRAM_OPTIONS_LIB) $(GSL_LIBS)

ngila_SOURCES = align.cpp matparser.cpp models.cpp ngila.cpp seqdb.cpp \
	seqparser.cpp


noinst_HEADERS = align.h matparser.h models.h ngila.cmds ngila.h seqdb.h \
	seqparser.h
