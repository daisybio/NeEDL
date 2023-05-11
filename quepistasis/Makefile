# find the flag using command 'python3.9-config --cflags --embed'
PYTHON_CFLAGS=-I/usr/include/python3.9 -I/usr/include/python3.9  -Wno-unused-result -Wsign-compare -g -fstack-protector-strong -Wformat -Werror=format-security  -DNDEBUG -g -fwrapv
# find the flags using command 'python3.9-config --ldflags --embed'
PYTHON_LDFLAGS=-L/usr/lib/python3.9/config-3.9-x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu -lpython3.9 -lcrypt -lpthread -ldl  -lutil -lm -lm

CPPFLAGS=-g -I header $(PYTHON_CFLAGS) 
LDFLAGS=-g $(PYTHON_LDFLAGS) 
LDLIBS=

test: build/cpu_sa.o build/cpu_sa_wrapper.o build/test.o build/snps_optimization.o build/python_wrapper.o
	g++ -o build/test build/cpu_sa_wrapper.o build/cpu_sa.o build/test.o build/snps_optimization.o build/python_wrapper.o $(LDFLAGS) $(LDLIBS)

build/python_wrapper.o: src/python_wrapper.cpp
	g++ $(CPPFLAGS) -c src/python_wrapper.cpp -o build/python_wrapper.o

build/cpu_sa.o: src/cpu_sa.cpp
	g++ $(CPPFLAGS) -c src/cpu_sa.cpp -o build/cpu_sa.o

build/cpu_sa_wrapper.o: src/cpu_sa_wrapper.cpp
	g++ $(CPPFLAGS) -c src/cpu_sa_wrapper.cpp -o build/cpu_sa_wrapper.o

build/snps_optimization.o: src/snps_optimization.cpp
	g++ $(CPPFLAGS) -c src/snps_optimization.cpp -o build/snps_optimization.o

build/main.o: src/main.cpp 
	g++ $(CPPFLAGS) -c src/main.cpp -o build/main.o

build/test.o: test/test.cpp 
	g++ $(CPPFLAGS) -c test/test.cpp -o build/test.o

clean:
	rm -f build/*
