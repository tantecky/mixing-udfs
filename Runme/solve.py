#!/usr/bin/env python3
import r2pipe

r2 = r2pipe.open('./runme', ['-d'])

while 'invalid' not in r2.cmd('s'):
    r2.cmd('ds;sr rip')

    json = r2.cmdj('pdj 1')

    if json is None:
        continue

    opcode = json[0]['opcode']

    if 'cmp byte [rcx]' in opcode:
        char = opcode.split(',')[1].strip()
        r2.cmd('wx {} @ rcx'.format(char))
        print(chr(int(char, 16)), end='', flush=True)

print()
