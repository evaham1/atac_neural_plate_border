#!/usr/bin/env python
 
print('Hello, World!')

import os
arr = os.listdir()
print(arr)

arr = os.listdir('./input/')
print(arr)

with open('readme.txt', 'w') as f:
    f.write('readme')