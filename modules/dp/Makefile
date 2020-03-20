all: dp_stdin

dp_stdin: dp_stdin_main.c dp.c dp.h
	$(CC) -O3 dp_stdin_main.c dp.c -o dp_stdin

clean:
	rm -f dp_stdin
