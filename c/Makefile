CC = cc
CFLAGS = -O3 -march=native -std=c99 -Wall

# Platform detection for framework linkage, FFTW3 and libyaml paths
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    FFTW_PREFIX := $(shell brew --prefix fftw 2>/dev/null || echo /opt/homebrew/opt/fftw)
    YAML_PREFIX := $(shell brew --prefix libyaml 2>/dev/null || echo /opt/homebrew/opt/libyaml)
    CFLAGS += -I$(FFTW_PREFIX)/include -I$(YAML_PREFIX)/include
    LDFLAGS = -lm -L$(FFTW_PREFIX)/lib -lfftw3 -L$(YAML_PREFIX)/lib -lyaml -framework Accelerate
else
    LDFLAGS = -lm -lfftw3 -lyaml -llapack -lblas
endif

OBJS = heat1dfun.o orbitfun.o fourier_solver.o yaml_config.o

all: heat1d

heat1d: heat1d.o $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

test_validate: test_validate.o heat1dfun.o orbitfun.o fourier_solver.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

test: test_validate
	./test_validate

clean:
	rm -f *.o heat1d test_validate

.PHONY: all test clean
