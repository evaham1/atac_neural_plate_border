#!/usr/bin/env python
 
print('Hello, World!')

import os
arr = os.listdir()
print(arr)

arr = os.listdir('./input/')
print(arr)

with open('readme.txt', 'w') as f:
    f.write('readme')

import matplotlib.pyplot as plt
plt.plot([0, 1, 2, 3, 4], [0, 3, 5, 9, 11])
plt.ylabel('Books Read')
plt.savefig('books_read.png')