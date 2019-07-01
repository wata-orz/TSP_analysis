#!/usr/bin/env python

k = 8

# Recursively enumerate all the connection patterns of k-moves and push them into patterns.
# M is the perfect matching on [2k], and M[i] is the mate of vertex i.
# -1 is a dummy value which means that the mate is not decided yet.
def enumerate_patterns(i, M, patterns):
    if i == 2 * k:
        if is_feasible(M):
            patterns.append(M[:])
    elif M[i] == -1:
        for j in range(i + 1, 2 * k):
            if M[j] == -1:
                (M[i], M[j]) = (j, i)
                enumerate_patterns(i + 1, M, patterns)
                M[i] = M[j] = -1
    else:
        enumerate_patterns(i + 1, M, patterns)

# Check the feasibility of the connection pattern M.
# Starting from vertex 0, we walk along the cycle and return whether the length is k.
def is_feasible(M):
    p = 0 # current vertex
    length = 0
    while p != 2 * k:
        length += 1
        p = M[p]
        if p % 2 == 0:
            p -= 1
        else:
            p += 1
    return length == k

# Decompose M into sequential swaps.
# The returned list [S_1, S_2, ..., S_c] is a partition of [k] such that each M[S_i] is a sequential swap.
# The list is sorted by |S_i|.
def sequential_swaps(M):
    seqs = []
    visited = [False] * (2 * k)
    for s in range(k):
        if not visited[s * 2]:
            # walk along the alternating cycle from 2s
            S = []
            p = s * 2
            while not visited[p]:
                S.append(p // 2)
                visited[p] = True
                p ^= 1
                visited[p] = True
                p = M[p]
            seqs.append(S)
    seqs.sort(key = len)
    return seqs

# Check whether M can be decomposed into two moves.
def is_reducible(M):
    seqs = sequential_swaps(M)
    c = len(seqs) # the number of sequential swaps
    for p in range(1, (1 << c) - 1): # try all the bipartitions of [c]
        M1 = M[:] # the swap consisting of sequential swaps not in p
        M2 = M[:] # the swap consisting of sequential swaps in p
        for i in range(c):
            if p >> i & 1:
                # delete the i-th sequential move from M1
                for j in seqs[i]:
                    (M1[j * 2], M1[j * 2 + 1]) = (j * 2 + 1, j * 2)
            else:
                # delete the i-th sequential move from M2
                for j in seqs[i]:
                    (M2[j * 2], M2[j * 2 + 1]) = (j * 2 + 1, j * 2)
        # If both the swaps are feasible, M is reducible.
        if is_feasible(M1) and is_feasible(M2):
            return True
    return False

# Swap i and j.
def swap(M, i, j):
    for p in range(2):
        mi = M[2 * i + p]
        mj = M[2 * j + p]
        (M[2 * i + p], M[mj]) = (mj, 2 * i + p)
        (M[2 * j + p], M[mi]) = (mi, 2 * j + p)


print("Enumerating connection patterns...")
patterns = []
enumerate_patterns(0, [-1] * (2 * k), patterns)
print("The number of connection patterns of {}-moves is {}.".format(k, len(patterns)))
patterns2 = []
for M in patterns:
    seqs = sequential_swaps(M)
    ks = list(map(len, seqs))
    if not ks == [2, 3, 3] or is_reducible(M):
        continue
    (i, j) = seqs[0]
    if i - 1 in seqs[1] and i + 1 in seqs[1] and j - 1 in seqs[2] and j + 1 in seqs[2]:
        patterns2.append(M)
    elif i - 1 in seqs[2] and i + 1 in seqs[2] and j - 1 in seqs[1] and j + 1 in seqs[1]:
        patterns2.append(M)
print("Among them, {} patterns meet the precondition in the lemma.".format(len(patterns2)))
for M in patterns2:
    seqs = sequential_swaps(M)
    ok = False
    for a in seqs[0]:
        if a - 1 in seqs[1]:
            A = 1
        else:
            A = 2
        if a - 2 in seqs[A] or a + 2 in seqs[A]:
            continue
        M1 = M[:]
        swap(M1, a, a + 1)
        M2 = M[:]
        swap(M2, a, a - 1)
        if (is_feasible(M1) or is_reducible(M1)) and (is_feasible(M2) or is_reducible(M2)):
            ok = True
    if not ok:
        print("counter example: {}".format(M))
        exit()
print("All of them satisfy the statement in the lemma.")
