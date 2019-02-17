OUT = main
OBJ_DIR = obj
DEP_DIR = dep
BOARD = %main_board.c
INTERREST := $(shell find . | grep \\.[hcS]$)
SRC = $(filter-out $(BOARD), $(filter %.c %.S, $(INTERREST)))
OBJ = $(addprefix $(OBJ_DIR)/, $(addsuffix .o, $(basename $(notdir $(SRC)))))
DEP = $(addprefix $(DEP_DIR)/, $(addsuffix .d, $(basename $(notdir $(INTERREST)))))
INC = $(addprefix -I , $(sort $(dir $(filter %.h, $(INTERREST)))))
VPATH = $(sort $(dir $(SRC)))

#DEBUG
#$(info $$SRC is [${SRC}])
#$(info $$OBJ is [${OBJ}])
#$(info $$DEP is [${DEP}])
#$(info $$INC is [${INC}])
#$(info $$VPATH is [${VPATH}])

# Konfiguration
CC = clang
CFLAGS = -Wall -Wextra -pedantic -O2 -std=c99
CPPFLAGS = $(INC)

# Regeln
.PHONY: all
all: main

-include $(DEP)

# generate dependencies
$(DEP_DIR)/%.d: %.c
	@set -e; rm -f $@; \
	$(CC) -MM $(CPPFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

# build objects
$(OBJ_DIR)/%.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

# build out
$(OUT): $(OBJ)
	$(CC) -o $@ $(OBJ) -lm

.PHONY: clean
clean:
	rm -f $(OUT)
	rm -f $(OBJ)
	rm -f $(DEP)
