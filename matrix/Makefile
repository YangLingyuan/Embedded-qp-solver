OUT = test
SRC = matrix.c test.c
OBJ = $(addsuffix .o, $(basename $(SRC)))

CC = clang
CFLAGS = -Wall -Wextra -pedantic -O2 -std=c99
MATH = -lm

$(OUT): $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $@ $(MATH)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -f $(OUT)
	rm -f $(OBJ)
