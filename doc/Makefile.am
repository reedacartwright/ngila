man_MANS = ngila.1

$(man_MANS): xcmd-doc

xcmd-doc: ngila.1pp

.1pp.1:
	$CPP $(CPPFLAGS) $(BOOST_CPPFLAGS) $(GSL_CFLAGS) -o $@ - < $<
