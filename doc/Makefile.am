man_MANS   = ngila.1
CLEANFILES = ngila.1
dist_man_MANS = mkman.pl

$(man_MANS): xcmd-doc

xcmd-doc: $(top_srcdir)/src/ngila.cmds
	$(srcdir)/mkman.pl < $(top_srcdir)/src/ngila.cmds > ngila.1
		

