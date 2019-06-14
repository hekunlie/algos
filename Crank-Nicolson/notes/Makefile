DATE:=    $(shell date +'%Y%m%d')
SRCS:=   notes.tex
default: notes

view: $(SRCS)
	latexmk -pdfxe -pvc $<

notes: $(SRCS)
	latexmk -pdfxe $<

clean:
	latexmk -C notes.tex
