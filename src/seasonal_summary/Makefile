#------------------------------------------------------------------------------
# Makefile for seasonal_summary executables
#------------------------------------------------------------------------------

# Set up compile options
CC = gcc
RM = rm -f
EXTRA = -Wall -g -O2

# Define the include files
INC1 = determine_max_extent.h
INC2 = generate_stack.h
INCDIR  = -I. -I$(XML2INC) -I$(ESPAINC)
NCFLAGS = $(EXTRA) $(INCDIR)

# Define the source code and object files
SRC1 = read_extent.c determine_max_extent.c
OBJ1 = $(SRC1:.c=.o)

SRC2 = read_xml.c generate_stack.c
OBJ2 = $(SRC2:.c=.o)

# Define the object libraries
LIB   = -L$(ESPALIB) -l_espa_raw_binary -l_espa_common -L$(XML2LIB) -lxml2 -lm

# Define the executable
EXE = determine_max_extent generate_stack

# Target for the executable
all: $(EXE)

determine_max_extent: $(OBJ1) $(INC1)
	$(CC) $(EXTRA) -o determine_max_extent $(OBJ1) $(LIB)

generate_stack: $(OBJ2) $(INC2)
	$(CC) $(EXTRA) -o generate_stack $(OBJ2) $(LIB)

install:
	install -d $(PREFIX)/bin
	install -m 755 $(EXE) $(PREFIX)/bin

clean:
	$(RM) *.o $(EXE)

$(OBJ1): $(INC1)
$(OBJ2): $(INC2)
.c.o:
	$(CC) $(NCFLAGS) -c $<

