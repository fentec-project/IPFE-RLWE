CC 		  = /usr/bin/gcc
CFLAGS  = -Wall -Wextra -O3 -fomit-frame-pointer -lcrypto -march=native 
CFLAGS2  = -Wall -Wextra -O3 -fomit-frame-pointer -march=native -std=c99 
CFLAGS3  = -Wall -Wextra -O3 -fomit-frame-pointer -march=native 
CLANG   = clang -march=native -O3 -fomit-frame-pointer -fwrapv -Qunused-arguments
RM 		  = /bin/rm


all: test/test_time \
     test/test_poly_mul \

SOURCES = crt.c modred.c poly_mul.c rlwe_mife.c rng.c sample.c
HEADERS = crt.h modred.h params.h poly_mul.h rlwe_mife.h rng.h sample.h


test/test_time: $(SOURCES) $(HEADERS) test/test_time.c
	$(CC) $(CFLAGS3) -o $@ $(SOURCES) test/test_time.c -lcrypto -lgmp

test/test_poly_mul: $(SOURCES) $(HEADERS) test/test_poly_mul.c
	$(CC) $(CFLAGS3) -o $@ $(SOURCES) test/test_poly_mul.c -lcrypto -lgmp

rng.o: rng.c
	$(CC) $(CFLAGS2) -c rng.c -lcrypto -o $@ 

.PHONY: clean test

test: 
	./test/test_time
	./test/test_poly_mul

clean:
	-$(RM) -r test/test_time
	-$(RM) -r test/test_poly_mul
