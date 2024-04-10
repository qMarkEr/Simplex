# Simplex method + sensitive analysis

## Basics
- **Simplex method** is a method for linear programming. It allows to solve optimization problems.
- **Sensitive analysis** determines the impact of changing the function shutdown coefficient (product prices in production planning problems) based on an automatic solution.

## Program 
- The `requirements.txt` file contains a description of the problem and values, that are parsed as input. 
- If you have a different model, change some values in `read_model.py` and the calls of 
    > equalizeMoreOrEqual(A, c, count)

    and
    > equalizeLessOrEquals(A, c, count)
 
    functions. 
  
    Also, use a 
  
    > solveX(A1, c1, m1, base_vectors)
  
   function in case you have imaginary variables, otherwise use 
   
   > SolveDefalt(A, c, basis, base_vectors)
   
    and pass the `determineBasis(A)` function call as a `basis` parameter.

## Usage
- Just run the:
    > python simplex.py