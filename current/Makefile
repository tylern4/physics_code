SUBDIRS = src cpp csv_maker golden_run
#python skim

all: clean $(SUBDIRS)

clean:
	-rm -f obj/*.o

$(SUBDIRS):
	$(MAKE) -C $@

.PHONY: clean $(SUBDIRS)
