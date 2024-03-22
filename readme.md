# code

We list the seven programs:

1_2.py: for an ode
2_3_Dirichlet.py: for a wave equation     
2_6_Neumann.py: for a wave equation 
2_6_Neumann_ghost.py: for a wave equation 
2_12_2D.py: for a 2-dimensional wave equation 
3_1.py: for a heat equation
burgers.py: for an one-dimensional burgers equation

# source

The first six programs are from the book:

@book{langtangen2017finite,
  title={Finite difference computing with PDEs: a modern software approach},
  author={Langtangen, Hans Petter and Linge, Svein},
  year={2017},
  publisher={Springer Nature}
}

The seventh program is from the scientist.

# explanation of program

For each program, `solver_scientist` includes the bug from the scientist.

`solver_mu1`  ~ `solver_mu10` are our ten mutants. 

The `main` function conduct our method. 