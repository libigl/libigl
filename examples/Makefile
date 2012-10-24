.PHONY: all
.NOTPARALLEL: all

DIRS=$(wildcard */)


all: $(DIRS) 
	for p in  $(DIRS); \
	do \
	echo "cd $$p" ; \
	$(MAKE) -C $$p; \
	done
	
clean: $(DIRS) 
	for p in  $(DIRS); \
	do \
	echo "cd $$p" ; \
	$(MAKE) -C $$p clean; \
	done
	#$(MAKE) -C $$p || exit $$?; \#$(MAKE) -C $$p || exit $$?; \
