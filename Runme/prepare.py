#!/usr/bin/env python3
import r2pipe

r2 = r2pipe.open('./runme', ['-w'])
end_enc = int(r2.cmd('S.~[0]'), 16)
#  r2.cmd('e scr.color=0')
r2.cmd('/c inc rcx')
res = r2.cmd('f*~hit[3]').splitlines()
enc_starts = [int(addr, 16) for addr in res]
keys = []

for enc_start in enc_starts:
    #  + 3 inc rcx
    enc = str(hex(enc_start + 3))
    key = r2.cmd('pd 1 @ ' + enc + '~cmp[5]')
    keys.append(key)

enc_starts.pop(0)
keys.pop()
#  print(enc_starts)
#  print(keys)

for enc_start, key in zip(enc_starts, keys):

    enc = str(hex(enc_start))
    num_bytes = str(hex(end_enc - int(enc, 16)))
    cmd = 'wox {} @ {}!{}'.format(key, enc, num_bytes)
    print(cmd)
    r2.cmd(cmd)
    #  break


