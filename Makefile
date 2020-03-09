
all: minimap2 dp_stdin

minimap2:
	cd modules/my_minimap2 && make && cd ../../

dp_stdin:
	cd modules/dp && make && cd ../../

clean:
	cd modules/my_minimap2 && make clean && cd ../../
	cd modules/dp && make clean && cd ../../
