CC ?= cc
CFLAGS ?= -O2
CFLAGS += -std=c11 -Wall -Wextra -Wpedantic
LDLIBS += -lm

.PHONY: all check clean

all: helios

helios: helios.c
	$(CC) $(CFLAGS) $< -o $@ $(LDLIBS)

check: helios
	./helios --help >/dev/null
	./helios 2026-08-17 | grep -q '06:39  Nascer do sol'
	./helios --local=35.67620,139.65030,Asia/Tokyo 2026-08-18 | grep -q 'Local temporário — 18/08/2026'
	! ./helios --local=91,0,Etc/UTC >/dev/null 2>&1
	! ./helios --ache a >/dev/null 2>&1
	./helios 2026-08-17 | grep -q '09:41  Nascer da lua'
	! ./helios 2026-02-30 >/dev/null 2>&1

clean:
	$(RM) helios
