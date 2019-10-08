CC =gcc
NO_ASS =-c
CFLAGS =-Wall
LINKS =-lm -ldislin -lgsl -lgslcblas
MAIN =main.c
MAIN_OUT =main.o
EXEC =main.out


all: compile tests clean

compile:
	$(info "COMPILING ...........")
	$(CC) $(CFLAGS) ${MAIN} ${NO_ASS} -o $(MAIN_OUT)
	$(CC) $(CFLAGS) $(MAIN_OUT) ${LINKS} -o $(EXEC)

tests:	
	$(info 'PROGRAM STARTING ...........')
	./$(EXEC)


clean:
	$(info 'PROGRAM DONE ...........')
	$(info 'CLEANING ...........')
	rm *.o*
