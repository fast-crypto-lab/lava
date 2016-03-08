CC=    gcc
CXX=   g++
#CXX=   ag++
LEX=   flex
YACC=  bison -d
LD=    g++

CFLAGS=  -O3 -std=gnu99 -msse2 -mpopcnt  -Werror
CXXFLAGS= -O3 -msse2 -mssse3 -mpopcnt   -pthread -fopenmp
INCPATH= -I/usr/local/include -I/opt/local/include -I/usr/include
LDFLAGS=  
LIBPATH= -L/usr/local/lib -L/opt/local/lib -L/usr/lib 
LIBS=    -lm -pthread -fopenmp -lm4ri


#TAR= 
EXE= f4_magma-test
OBJ= f4_algo.o 
CSRC= $(wildcard *.cpp)


ifdef DEBUG
        CFLAGS+=  -DDEBUG_MODE
        CXXFLAGS+= -DDEBUG_MODE
endif

ifdef GPROF
	CFLAGS += -pg
	CXXFLAGS += -pg
	LDFLAGS += -pg
endif

all: $(EXE)


$(EXE): $(OBJ) $(EXE).o

%-test: $(OBJ) %-test.o
	$(LD) $(LDFLAGS) $(LIBPATH) -o $@ $^ $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) $(INCPATH) -c $<

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCPATH) -c $<



clean:
	rm *-test *.o 
