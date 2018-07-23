SUBDIRS = src cpp csv_maker golden_run skim
#python skim

all: $(SUBDIRS)

clean:
	-rm -f obj/*.o

$(SUBDIRS):
	$(MAKE) -C $@

.PHONY: clean $(SUBDIRS)
