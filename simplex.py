import scipy as sc
import read_model as rm

def nice_print(A, c, basis, deltas):
    print('-'*(len(A[0] * 8) + 6))
    print("B", end='\t')
    print("CB", end='\t')
    for i in range(len(A[0])):
        print("P", i, sep = "", end = '\t')
    print()
    for i in range(len(A)):
        print(basis[i], end='\t')
        print(c[basis[i]], end='\t')
        for j in A[i]:
            print(round(j, 2), end='\t')
        print()
    print(end='\t\t')
    for i in deltas:
        print(round(i, 2), end='\t')
    print()

def nice_print_X(A, c, basis, deltas):
    print('-'*(len(A[0] * 8) + 6))
    print("B", end='\t')
    print("CB", end='\t')
    for i in range(len(A[0])):
        print("P", i, sep = "", end = '\t')
    print()
    for i in range(len(A)):
        print(basis[i], end='\t')
        print(c[basis[i]], end='\t')
        for j in A[i]:
            print(round(j, 2), end='\t')
        print()
    print(end='\t\t')
    for i in deltas:
        print(round(i[1], 2), end='\t')
    print()
    print(end='\t\t')
    for i in deltas:
        print(round(i[0], 2), end='\t')
    print()

def isBasis(vector):
    return (vector.count(1) == 1 and vector.count(0) == len(vector) - 1) or (vector.count(1.0) == 1 and vector.count(0.0) == len(vector) - 1)

def addVector(vector, A, c):
    for (i, j) in zip(A, vector):
        i.append(j)
    c.append(0)

def addVectorOmega(vector, A, c):
    for (i, j) in zip(A, vector):
        i.append(j)
    c.append(1)

def getVector(A, index):
    return [i[index] for i in A]

def equalizeLessOrEqual(A, c, count):
    for i in range(count):
        new_basis = [0 for _ in range(len(A))]
        new_basis[i] = 1
        addVector(new_basis, A, c)

def equalizeMoreOrEqual(A, c, count):
    for i in range(count):
        new_basis = [0 for _ in range(len(A))]
        new_basis[len(new_basis) - count + i] = -1
        addVector(new_basis, A, c)
    m_ind = []
    for i in range(count):
        new_basis = [0 for _ in range(len(A))]
        new_basis[len(new_basis) - count + i] = 1
        addVectorOmega(new_basis, A, c)
        m_ind.append(len(A[0]))
    return m_ind

def dot(x, y):
    return sum([i * j for (i, j) in zip(x, y)])

def CalcDeltas(A, c, basis):
    d = [0 for _ in range(len(A[0]))]
    CB = []
    for j in basis:
        CB.append(c[j])

    for j in range(len(d)):
        temp_vector = getVector(A, j)
        d[j] = dot(temp_vector, CB) - c[j]
    return d

def replace(A, basis, delta):
    src_index = delta.index(max(delta[1::]))
    vec = getVector(A, src_index)
    dst_index = 0
    ratio = 999999
    for i in range(len(vec)):
        if vec[i] > 0 and A[i][0] > 0:
            if A[i][0] / vec[i] < ratio:
                ratio = A[i][0] / vec[i]
                dst_index = i
    basis[dst_index] = src_index
    return src_index

def substract(A, i, j):
    divisor = A[j][i]
    A[j] = [i / divisor for i in A[j]]
    for index in range(len(A)):
        if index != j:
            multiplier = A[index][i] / A[j][i]
            A[index] = [elem2 + elem1 * -multiplier for (elem1, elem2) in zip(A[j], A[index])]

def determineBasis(A):
    basis = []
    for i in range(1, len(A[0])):
        vec = getVector(A, i)
        if isBasis(vec):
            basis.append(i)
    return basis

def CalcDoubleDeltas(A, c, basis, m_indices):
    d = [[0, 0] for _ in range(len(A[0]))] # d[1] * omega + d[0]
    CB = []
    for j in basis:
        CB.append(c[j])

    for j in range(len(d)):
        temp_vector = getVector(A, j)
        for index in m_indices:
            temp_vector[basis.index(index)] = 0
        d[j][0] = dot(temp_vector, CB) - (c[j] if j not in m_indices else 0)

    for j in range(len(d)):
        temp_vector = getVector(A, j)
        for i in range(len(temp_vector)):
            if basis[i] not in m_indices:
                temp_vector[i] = 0
        d[j][1] = dot(temp_vector, CB) - (c[j] if j in m_indices else 0)
    
    return d

def replaceOmega(A, c, basis, m_indices, delta):
    src_index = 0
    max_temp = 0
    for i in range(1, len(delta)):
        if sum(delta[i]) > max_temp and delta[i][1] != 0:
            src_index = i
            max_temp = sum(delta[i])
    vec = getVector(A, src_index)
    dst_index = -1
    ratio = 999999
    for i in range(len(vec)):
        if basis[i] in m_indices:
            if vec[i] > 0 and A[i][0] > 0:
                if A[i][0] / vec[i] < ratio:
                    ratio = A[i][0] / vec[i]
                    dst_index = i
    
    c[basis[dst_index]] = 0
    for i in range(len(A)):
        A[i][basis[dst_index]] = 0
    
    m_indices.remove(basis[dst_index])
    basis[dst_index] = src_index
    return src_index

def SolveDefault(A, c, basis, base_vectors):
    deltas = CalcDeltas(A, c, basis)
    nice_print(A, c, basis, deltas)

    # main loop
    while any(x > 0 for x in deltas[1:]):
        new_vec = replace(A, basis, deltas)
        substract(A, new_vec, basis.index(new_vec))
        deltas = CalcDeltas(A, c, basis)
        nice_print(A, c, basis, deltas)

    print("Optimum:", -deltas[0])
    print("Result vector: {", end=' ')
    for i in range(1, base_vectors + 1):
        if i in basis:
            print(A[basis.index(i)][0], end=' ')
        else:
            print(0, end=' ')
    print("}")

    
def solveX(A1, c1, m1, base_vectors):
    try:
        basis = determineBasis(A1)
        deltas = CalcDoubleDeltas(A1, c1, basis, m1)
        nice_print_X(A1, c1, basis, deltas)
        while len(m1) != 0:
            new = replaceOmega(A1, c1, basis, m1, deltas)
            substract(A1, new, basis.index(new))
            deltas = CalcDoubleDeltas(A1, c1, basis, m1)
            nice_print_X(A1, c1, basis, deltas)
        SolveDefault(A1, c1, basis, base_vectors)
        for i in range(3):
            temp = input()
            if temp != 'q':
                c1[i] = int(temp)
            else:
                return
        SolveDefault(A1, c1, basis, base_vectors)
    except:
        print("smth wrong with your model")

# input data
A_input, b_input, c_input = rm.fillMatrices()

# check
# res = sc.optimize.linprog(A_ub= [[2, 5, 1], [5, 30, 2]], b_ub=[250, 1000], c=[-105, -450, -13], bounds=(10, None), method='simplex')
# print(res)

# max index of default vector
base_vectors = len(A_input[0])

# add fake vectors
equalizeLessOrEqual(A_input, c_input, 2)
m_input = equalizeMoreOrEqual(A_input, c_input, 3)

#add P0
for (i, j) in zip(A_input, b_input):
    i.insert(0, j)
c_input.insert(0, 0)

# function-driver
solveX(A_input, c_input, m_input, base_vectors)