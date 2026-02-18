CC = cc
CFLAGS = -O3 -march=native -std=c99 -Wall

# Platform detection for framework linkage and FFTW3 paths
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    FFTW_PREFIX := $(shell brew --prefix fftw 2>/dev/null || echo /opt/homebrew/opt/fftw)
    CFLAGS += -I$(FFTW_PREFIX)/include
    LDFLAGS = -lm -L$(FFTW_PREFIX)/lib -lfftw3 -framework Accelerate
else
    LDFLAGS = -lm -lfftw3 -llapack -lblas
endif

OBJS = heat1dfun.o orbitfun.o fourier_solver.o

all: heat1d_moon

heat1d_moon: heat1d_moon.o $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

test_validate: test_validate.o $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

test: test_validate
	./test_validate

clean:
	rm -f *.o heat1d_moon test_validate

.PHONY: all test clean
