#!/usr/bin/env python3
import r2pipe

r2 = r2pipe.open('./runme', ['-d'])
key = []

i = 0

while True:
    r2.cmd('ds;sr rip')
    line = r2.cmd('pd 1')
    if 'cmp byte' in line:
        print(line)
        letter = r2.cmd('pd 1~cmp[5]')
        print(letter)
        r2.cmd('wx {} @ rcx'.format(letter))
        print(r2.cmd('p8 1 @ rcx'))
        i += 1



