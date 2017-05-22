BITS 64

section .mytext progbits alloc exec write align=16
global _start:

%macro sys_exit 1
  mov eax, 0x3c
  mov edi, %1
  syscall
%endmacro

%macro encrypt_below 1
  ; key
  mov eax, %1
  ; length of a block to encrypt
  mov edi, enc_end - %%enc_start
  xor esi, esi
%%do_for:
  cmp esi, edi
  je %%enc_start
  xor [%%enc_start + esi], al
  inc esi
  jmp %%do_for
%%enc_start:
%endmacro

_start:
; :> wox 0x58 @0x006000cb!0xc
  encrypt_below 'X'
  sys_exit 0
enc_end:
  sys_exit 0
