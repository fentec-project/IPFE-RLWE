CC 		  = /usr/bin/gcc
CFLAGS  = -Wall -Wextra -O3 -fomit-frame-pointer -march=native -mavx2
RM 		  = /bin/rm

all: test/test_time

SOURCES = arith_rns.c crt.c randombytes.c gauss.c ntt.c rlwe_sife.c sample.c aes256ctr.c
HEADERS = arith_rns.h crt.h randombytes.h gauss.h ntt.h rlwe_sife.h sample.h aes256ctr.h consts.h params.h


test/%: test/%.c $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) -o $@ $(SOURCES) $< -lgmp

.PHONY: clean test

test: 
	./test/test_time

clean:
	-$(RM) -rf test/test_time
