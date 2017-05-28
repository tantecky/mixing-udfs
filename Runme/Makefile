AS := nasm
ASM_FILES := runme.asm
OBJ_FILES := $(ASM_FILES:.asm=.o)
BIN_FILE := runme

$(BIN_FILE): $(OBJ_FILES)
	$(LD) -s -z noexecstack -o $(BIN_FILE) $(OBJ_FILES)
	python3 ./prepare.py

$(OBJ_FILES): $(ASM_FILES)
	$(AS) -f elf64 -o $@ $<

