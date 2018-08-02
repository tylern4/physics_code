SUBDIRS = cpp csv_maker golden_run skim
#python skim

all: lib $(SUBDIRS)

lib:
	$(MAKE) -j8 -C src

clean:
	-rm -f obj/*.o

$(SUBDIRS):
	$(MAKE) -C $@

.PHONY: clean $(SUBDIRS)
