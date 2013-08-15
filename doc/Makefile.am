
LATEX=@LATEX@
DVIPS=@DVIPS@
DVIPDF=@DVIPDF@
DOXYGEN=@DOXYGEN@

SUFFIXES = .tex .dvi .ps .pdf

.tex.dvi:
	$(LATEX) $<
	$(LATEX) $<

.dvi.ps:
	$(DVIPS) $< -o $@

.dvi.pdf:
	$(DVIPDF) $< 

texfile = UserGuide.tex

if LATEXTEST
dvifile = UserGuide.dvi
else
dvifile = 
endif

if DVIPSTEST
psfile = UserGuide.ps
else
psfile = 
endif

if DVIPDFTEST
pdffile = UserGuide.pdf
else
pdffile = 
endif


dist_doc_DATA_files = $(dvifile) $(psfile) $(pdffile) 

doc_DATA = $(dist_doc_DATA_files)

if DOXYGENTEST
doxytest = doxygen
else
doxytest = 
endif

doc_DATA += $(doxytest)

docdir=@datadir@/docs
mandir=@mandir@

doxygen:
	$(DOXYGEN)

install-docDATA:
	mkdir -p $(docdir)
	test `echo $(doc_DATA_files) | awk '{print NF}'` -ne 0 && cp $(doc_DATA_files) $(docdir) || echo
	test -d ./html && rm -rf crpropa_html &&  mv ./html crpropa_html || echo; 
	rm -rf $(docdir)/crpropa_html
	test -d ./crpropa_html && cp -r ./crpropa_html $(docdir)/ || echo
	test -d $(mandir) || mkdir $(mandir)
	test -d $(mandir)/man3 || mkdir $(mandir)/man3
	cp man/man3/* $(mandir)/man3/
	
uninstall-docDATA:
	cd $(docdir) ; rm -rf DraftUserGuide.ps *.pdf *.dvi crpropa_html/ 
	
clean:
	rm -rf html crpropa_html
	rm -rf *.toc *.log *.lof *.aux UserGuide.ps *.dvi *.pdf

distclean: clean uninstall-docDATA

