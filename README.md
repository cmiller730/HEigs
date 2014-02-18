HEigs
=====

Haskell interface to arpack for large sparse eigenvalue problems.

ARPACK is a Fortran code for computing a few eigenpairs associated with large sparse linear systems.
This package wraps a subset of ARPACK's functionality and attempts to deliver something similar to
scipy or MATLAB's eigs functions. 

To solve an eigen system Ax = \lambda x the user needs to define an ArpackLinearOp

type ArpackLinearOp = (SV.IOVector CDouble -> SV.IOVector CDouble -> IO ())

This operator should overwrite the second vector with the matrix times the first vector. To compute 
the eigenvalues call

eigs :: ArpackLinearOp -> ProblemDim -> Which -> NumEV -> Tolerance -> 
        MaxIter -> IO (Bool, [(Complex Double, V.Vector (Complex Double))])
        
Where, 

type ProblemDim = Int         -- The size of the linear system.
data Which = LM | SM deriving Show  -- Which eigenvalues to compute Largest magnitude or Smallest magnitude
                                    -- will support additional modes in future.
type NumEV = Int              -- The number of eigenpairs to compute.
type Tolerance = Double       -- The error tolerance. Setting to 0 uses the machine eps.
type MaxIter = Int            -- The maximum number of Arnoldi iterations.

The return contains (False, []) if ARPACK was unsuccessful.
Otherwise the second element in the tuple contains a list of the computed NumEV eigenvalues and eigenvectors.

History
=======

0.0.1 First working interface to ARPACK. Can solve Ax = lambda x for largest and smallest magnitude 
      eigenproblems.
      
Future Plans
============

Support additional values for Which.
Write a second wrapper for ARPACK's symmetric driver.

I usually only write what I need for today. If you need to access other portions of ARPACKs functionality
please let me know and I'll do my best to support it.

Chris Miller
cmiller730@gmail.com
