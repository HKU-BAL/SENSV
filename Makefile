
all: minimap2 dp_stdin

minimap2:
	cd my_minimap2 && make && cd ..

dp_stdin:
	cd dp && make && cd ..

clean:
	cd my_minimap2 && make clean && cd ..
	cd dp && make clean && cd ..
