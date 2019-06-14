CC         =  gcc #mpicc
OPTS       =  
LIBS       =  

EXEC       =  test
SRCS       =  $(wildcard *.c)
MY_INCL    =  $(wildcard *.h)
OBJS       =  $(SRCS:.c=.o)

$(EXEC): $(OBJS)
	$(CC) $(OPTS) $(OBJS) $(LIBS) -lm -o $(EXEC)

%.o:%.c $(MY_INCL) Makefile
	$(CC) $(OPTS) $(DEBUG) $(INCL) -c $< -o $@

.PHONY: clean

clean:
	-rm  $(EXEC)  $(OBJS)
