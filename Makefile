CC = g++
CFLAGS = -std=c++0x -lMinuit
LFLAGS = -I. -lm
LIBSMAIN = `root-config --cflags --evelibs`
LIBS = `root-config  --cflags --evelibs`

SOURCES  := $(wildcard *.C)
INCLUDES := $(wildcard *.h)
OBJECTS  := $(SOURCES:%.C=%.o)

all: example

example: $(OBJECTS)
	@$(CC) $(OBJECTS) $(LFLAGS) -o example.x $(LIBSMAIN)

$(OBJECTS): %.o : %.C
	@$(CC) $(CFLAGS) -c $< -o $@ $(LIBSMAIN)


clean:
	@rm -f $(OBJECTS)




