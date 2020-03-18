import numpy as np
#function HammingG
#input: a number r
#output: G, the generator matrix of the (2^r-1,2^r-r-1) Hamming code
def hammingGeneratorMatrix(r):
    n = 2**r-1
    
    #construct permutation pi
    pi = []
    for i in range(r):
        pi.append(2**(r-i-1))
    for j in range(1,r):
        for k in range(2**j+1,2**(j+1)):
            pi.append(k)

    #construct rho = pi^(-1)
    rho = []
    for i in range(n):
        rho.append(pi.index(i+1))

    #construct H'
    H = []
    for i in range(r,n):
        H.append(decimalToVector(pi[i],r))

    #construct G'
    GG = [list(i) for i in zip(*H)]
    for i in range(n-r):
        GG.append(decimalToVector(2**(n-r-i-1),n-r))

    #apply rho to get Gtranpose
    G = []
    for i in range(n):
        G.append(GG[rho[i]])

    #transpose    
    G = [list(i) for i in zip(*G)]

    return G


#function decimalToVector
#input: numbers n and r (0 <= n<2**r)
#output: a string v of r bits representing n
def decimalToVector(n,r): 
    v = []
    for s in range(r):
        v.insert(0,n%2)
        n //= 2
    return v

# Problem 1
def message(a):
    l = len(a)
    r = 2
    m =[]

    while 2**r - 2*r - 1 < l:
        r += 1

    k = 2**r - r - 1    

    l_binary = decimalToVector(l,r)

    for bit in l_binary:
        m.append(int(bit))

    for bit in a:
        m.append(int(bit))

    while len(m) < k:
        m.append(0)

    return m
 

# Problem 2
def hammingEncoder(m):
    l = len(m)
    r = 2

    while True :
        x = 2**r - r - 1
        if x > l:
            return []
        if x < l:
            r += 1
        if x == l: 
            break

    G_matrix = np.matrix(hammingGeneratorMatrix(r))
    m_matrix = np.matrix(m)

    output = m_matrix*G_matrix

    output = output.tolist()[0]
    final = []

    for bit in output:
        bit = bit%2
        final.append(bit)
    
    return final    

#Problem 3
def hammingDecoder(v):
    l = len(v)
    r = 2

    while True :
        x = 2**r - 1
        if x > l:
            return []
        if x < l:
            r += 1
        if x == l: 
            break

    H_matrix_transpose = []

    if r==2:
        j=3
    else:
        j=r

    for i in range(1,x+1):
        empty = []
        bit_rep = decimalToVector(i,j)
        for bit in bit_rep:
            empty.append(bit)
        H_matrix_transpose.append(empty)  

    v_ = v

    syndrome = (np.matrix(v_) * np.matrix(H_matrix_transpose)).tolist()[0]


    final = []
    for bit in syndrome:            
        bit = bit%2
        final.append(bit)
        
    check_zero = all(elem == 0 for elem in final)
    if check_zero:
        return v
   
    final_string = ''.join(str(bit) for bit in final)
    final = int(final_string, 2)


    if v[final-1] == 1:
        v[final-1] = 0

    elif v[final-1] == 0:
        v[final-1] = 1

    return v 


#Problem 4
def messageFromCodeword(c):

    l = len(c)
    r = 2
    output = []

    while True :
        x = 2**r - 1
        if x > l:
            return []
        if x < l:
            r += 1
        if x == l: 
            break

    positions = []

    for i in range (0,r):
        positions.append(2**i)

    i = 1 
    for element in c:
        if i not in positions:
            output.append(c[i-1])
        i +=1

    return output

#Problem 5
def dataFromMessage(m):
    l = len(m)
    r = 2
    while True :
        x = 2**r -r - 1
        if x > l:
            return []
        if x < l:
            r += 1
        if x == l: 
            break
    k = 2**r - r - 1 
 

    l_raw = []
    for element in m[0:r]:
        l_raw.append(element)
    out = 0
    for bit in l_raw:
            out = (out << 1) | bit
    l_raw = out
   
    if l_raw > k-r:
        return []

    raw_data = []
    for element in m[r:l_raw+r]:
        raw_data.append(element)

    return raw_data

#Problem 6
def repetitionEncoder(m,n):

    output = []
    m_int = m[0]
    
    for i in range(0,n):
        output.append(m_int)

    return output     


#Problem 7
def repetitionDecoder(v):
    count_0 = 0
    count_1 = 0
    for entry in v:
        if entry == 0:
            count_0 +=1
        if entry == 1:
            count_1 +=1
    if count_0 == count_1:
        return []
    elif count_0 > count_1:
        return [0]
    else:
        return [1] 
