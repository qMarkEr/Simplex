from simplex import *

# input data
try:
    A_input, b_input, c_input, more, less = rm.fillMatrices()
    base_vectors = len(A_input[0])

    # add fake vectors
    equalizeLessOrEqual(A_input, c_input, less)
    m_input = equalizeMoreOrEqual(A_input, c_input, more)

    #add P0
    for (i, j) in zip(A_input, b_input):
        i.insert(0, j)
    c_input.insert(0, 0)

    # function-driver
    solveX(A_input, c_input, m_input, base_vectors)
except:
    print(":(")