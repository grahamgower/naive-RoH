TARGET=hom_windows
HTS=../htslib
CFLAGS=-Wall -O2 -g -I$(HTS)
#LDLIBS=-L$(HTS) -Wl,-static -lhts -lz -Wl,-Bdynamic -pthread
LDLIBS=-L$(HTS) -lhts
CC=gcc

$(TARGET): $(TARGET).o
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)

clean:
	rm -f $(TARGET) $(TARGET).o
