BITS 64
%assign	STDIN	0
%assign	STDOUT 1

section .rodata
msg_input db '> ',0
msg_input_len equ $-msg_input
msg_err db `Wrong!\n`,0
msg_err_len equ $-msg_err


section .mytext progbits alloc exec write align=16
global _start:

%macro sys_exit 1
  mov eax, 0x3c
  mov edi, %1
  syscall
%endmacro

%macro sys_write 2
  mov eax, 0x1
  ; stdout
  mov edi, STDOUT
  ; buffer
  mov rsi, %1
  ; len
  mov edx, %2
  syscall
%endmacro

%macro sys_read 2
  mov eax, 0x0
  ; stdout
  mov edi, STDIN
  ; buffer
  mov rsi, %1
  ; len
  mov edx, %2
  syscall
%endmacro

%macro encrypt_below 0
  ; key
  mov rax, [rcx]
  ; length of a block to encrypt
  mov edi, enc_end - %%enc_start
  xor esi, esi
%%do_for:
  cmp esi, edi
  jz %%enc_start
  xor [%%enc_start + esi], al
  inc esi
  jmp %%do_for
%%enc_start:
%endmacro

%macro test_char 1
  ; provided input/key
  cmp [rcx], byte %1
  jz %%good
  sys_write msg_err, msg_err_len
  sys_exit 1
%%good:
  inc rcx
%endmacro

_start:
; :> wox 0x58 @0x006000cb!0xc
  sys_write msg_input, msg_input_len
  sub rsp, 16
  sys_read rsp, 16
  ; pointer to input
  mov rcx, rsp
  test_char 'c'
  encrypt_below
  test_char 'r'
  encrypt_below
  test_char 'y'
enc_end:
  sys_exit 0

