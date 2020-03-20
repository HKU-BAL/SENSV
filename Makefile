
all: minimap2 dp_stdin dpmodule

minimap2:
	cd modules/my_minimap2 && make && cd ../../

dp_stdin:
	cd modules/dp && make && cd ../../

dpmodule:
	cd modules/dpmodule && $(PYTHON) setup.py install && cd ../../

clean:
	cd modules/my_minimap2 && make clean && cd ../../
	cd modules/dp && make clean && cd ../../
	cd modules/dpmodule && rm -r build && cd ../../
