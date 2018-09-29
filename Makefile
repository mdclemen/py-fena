
funcs.so: funcs.o
	$(CC) -fPIC funcs.o -shared -o funcs.so

funcs.o:
	$(CC) -c -fPIC funcs.c

clean:
	rm *.o *.so *~

all: funcs.so
