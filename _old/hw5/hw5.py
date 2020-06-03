#!/usr/bin/env python
# coding: utf-8

# In[9]:


get_ipython().system('jupyter nbconvert --to script hw5.ipynb')
get_ipython().system('jupyter nbconvert --to pdf hw5.ipynb')


# In[13]:


from typing import *
import numpy as np
import bitstring as bs


# In[14]:


# max_len determines the maximum length in the second dimension of Python list
def max_len(items: List) -> int:
    maxLen = 0
    for item in items:
        if (itemLen := len(item)) > maxLen:
            maxLen = itemLen
    return maxLen

# uniform_octree_len expands octree levels to where the depths are uniform
def uniform_octree_len(rSequences: List, maxLen: int) -> List:
    shortLenMaskFunc = lambda rSeq: len(rSeq) != maxLen
    while any(shortLenMask := list(map(shortLenMaskFunc, rSequences))):
        i = shortLenMask.index(True)
        rSequences.extend([rSequences[i].copy() + [str(digit)] for digit in range(8)])
        del rSequences[i]
    
    return sorted(rSequences)

# lookup_table for Table 2 from Ahuja and Nash
def lookup_table(label: Union[str, int], direction: Union[str, int]) -> str:
    if direction not in {"x", "y", "z", "0", "1", "2", 0, 1, 2}:
        raise IndexError("direction must be in {'x', 'y', 'z', '0', '1', '2', 0, 1, 2}")
    if type(label) == str:
        label = int(label)
    if label > 7 or label < 0:
        raise IndexError("label must be between 0 and 7")
    
    if direction in ("x", "y", "z"):
        direction = ("x", "y", "z").index(direction)
    elif direction in ("0", "1", "2"):
        direction = int(direction)
    
    table = (
        (1, 2, 4), 
        (10, 3, 5),  
        (3, 10, 6),  
        (12, 11, 7),  
        (5, 6, 10),  
        (14, 7, 11), 
        (7, 14, 12), 
        (16, 15, 13)
    )
    return str(table[label][direction])

# octree_displacement returns the rSequence after it has been translated in `direction` by the `displacement` (binary/2's comp) using the `lookup_table`
def octree_displacement(rSeq: str, displacement: str, direction: Union[int, str]) -> List:
    # copy so we don't mutate outside lists
    displacement = displacement.copy() 
    rSeq = rSeq.copy()
    
    for i in reversed(range(len(displacement))):
        if int(displacement[i]) == 2:
            displacement[i] = str(0)
            if i != 0:
                displacement[i-1] = str(int(displacement[i-1]) + 1)
        
        if int(displacement[i]) == 0:
            continue
        
        tableVal = lookup_table(rSeq[i], direction)
        if int(tableVal) < 10:
            rSeq[i] = tableVal
        else:
            rSeq[i] = tableVal[-1]
            if i != 0:
                displacement[i-1] = str(int(displacement[i-1]) + 1)
    
    return rSeq


# In[15]:


# read in A and B from the text files (copied from .pdf document)
with open("ObjectA.txt", "rt") as f:
    A  = [rSeq.strip() for rSeq in f.read().split(",")]

with open("ObjectB.txt", "rt") as f:
    B  = [rSeq.strip() for rSeq in f.read().split(",")]

# split into letter lists & drop the 'r'
A = [list(rSeq)[1:] for rSeq in A]
B = [list(rSeq)[1:] for rSeq in B]

# find the length needed to expand to
if (ALen := max_len(A)) >= (BLen := max_len(B)):
    maxLen = ALen
elif (ALen := max_len(A)) < (BLen := max_len(B)):
    maxLen = BLen

# make uniform lengths
A = uniform_octree_len(A, maxLen)
B = uniform_octree_len(B, maxLen)

# x and y binary/2's comp representation
x = list(str(bs.Bits(int=-5, length=maxLen).bin))
y = list(str(bs.Bits(int=48, length=maxLen).bin))

# translate B
BPrime = []
for b in B:
    BPrimeX = octree_displacement(b, x, "x")
    BPrimeXY = octree_displacement(BPrimeX, y, "y")
    BPrime.append(BPrimeXY)

# test two results
assert("".join(BPrime[7]) == "52220026")
assert("".join(BPrime[-1]) == "70007036")

# write it out to file
print(f"B': ")
with open("BPrime.txt", "wt") as f:
    for b in BPrime:
        out = f"r{''.join(b)}"
        print(out)
        f.write(f"{out}\n")
    
# need to convert A and BPrime to hashable containers
A = set([tuple(a) for a in A])
BPrime = set([tuple(b) for b in BPrime])

# find intersection and volume
intersection = A & BPrime
volSmallestOctree = (100 / (2**maxLen) ** 3)
volTotal = len(intersection) * volSmallestOctree
print(f"Intersection volume: {volTotal} cm^3")


# In[ ]:




