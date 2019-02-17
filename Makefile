OUT = main
OBJ_DIR = obj
DEP_DIR = dep
INTERREST := $(shell find . | grep \\.[hcS]$)
SRC = $(filter-out $(BOARD), $(filter %.c %.S, $(INTERREST)))
OBJ = $(addprefix $(OBJ_DIR)/, $(addsuffix .o, $(basename $(notdir $(SRC)))))
INC = $(addprefix -I , $(sort $(dir $(filter %.h, $(INTERREST)))))
VPATH = $(sort $(dir $(SRC)))

#DEBUG
#$(info $$SRC is [${SRC}])
#$(info $$OBJ is [${OBJ}])
#$(info $$INC is [${INC}])
#$(info $$VPATH is [${VPATH}])

# Konfiguration
CC = clang
DEFS = -DNO_BOARD
CFLAGS = -Wall -Wextra -pedantic -O2 -std=c99 $(DEFS)
CPPFLAGS = $(INC)

# Regeln
.PHONY: all
all: main

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
