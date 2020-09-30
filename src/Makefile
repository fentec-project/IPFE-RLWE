CC 		  = /usr/bin/gcc
CFLAGS  = -Wall -Wextra -O3 -fomit-frame-pointer -march=native -mavx2
RM 		  = /bin/rm

all: test/test_time

SOURCES = crt.c randombytes.c modred.c gauss.c poly_mul.c mod_red.c ntt.c rlwe_mife.c sample.c aes256ctr.c
HEADERS = crt.h randombytes.h modred.h gauss.h consts.h mod_red.h ntt.h params.h poly_mul.h rlwe_mife.h sample.h aes256ctr.h


test/test_time: $(SOURCES) $(HEADERS) test/test_time.c
	$(CC) $(CFLAGS) -o $@ $(SOURCES) test/test_time.c -lgmp

.PHONY: clean test

test: 
	./test/test_time

clean:
	-$(RM) -r test/test_time
