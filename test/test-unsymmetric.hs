import Test.Tasty
import Test.Tasty.SmallCheck as SC
import Test.Tasty.QuickCheck as QC
import Test.Tasty.HUnit

import qualified Data.Vector.Storable.Mutable as MV
import qualified Data.Vector.Storable as V
import Numeric.LinearAlgebra.UnsymmetricEigen
import Control.Monad
import Foreign.C.Types
import System.IO.Unsafe
import Data.Complex

import Data.List
import Data.Ord

main = defaultMain tests

tests :: TestTree
tests = testGroup "Tests" [unitTests]

{-
properties :: TestTree
properties = testGroup "Properties" [scProps, qcProps]

scProps = testGroup "(checked by SmallCheck)"
  [ SC.testProperty "sort == sort . reverse" $
      \list -> sort (list :: [Int]) == sort (reverse list)
  , SC.testProperty "Fermat's little theorem" $
      \x -> ((x :: Integer)^7 - x) `mod` 7 == 0
  -- the following property does not hold
  , SC.testProperty "Fermat's last theorem" $
      \x y z n ->
        (n :: Integer) >= 3 SC.==> x^n + y^n /= (z^n :: Integer)
  ]

qcProps = testGroup "(checked by QuickCheck)"
  [ QC.testProperty "sort == sort . reverse" $
      \list -> sort (list :: [Int]) == sort (reverse list)
  , QC.testProperty "Fermat's little theorem" $
      \x -> ((x :: Integer)^7 - x) `mod` 7 == 0
  -- the following property does not hold
  , QC.testProperty "Fermat's last theorem" $
      \x y z n ->
        (n :: Integer) >= 3 QC.==> x^n + y^n /= (z^n :: Integer)
  ]
-}

unitTests = testGroup "Unit tests"
  [ testCase "Diagonal Matrix LM" $ diagMatLMCheck @?= True,
    testCase "Diagonal Matrix SM" $ diagMatSMCheck @?= True,
    testCase "Non-Symmetric Sanity" $ complexEigsCheck @?= True
  ]
  
diagMat :: ArpackLinearOp
diagMat x y = do
  MV.copy y x
  forM [0..999] (\n -> do x' <- MV.read x n
                          let lambda = fromIntegral (n+1) :: CDouble
                              y' = x'*(fromIntegral (n+1))
                          MV.write y n y'
                )
  return ()

complexEigs :: ArpackLinearOp
complexEigs x y = do
  x' <- V.freeze x
  MV.write y 0 (3*((V.!) x' 0) + (-2)*((V.!) x' 1))
  MV.write y 1 (4*((V.!) x' 0) + (-1)*((V.!) x' 1))
  MV.write y 2 (2*((V.!) x' 2))
  MV.write y 3 (1*((V.!) x' 3))
  MV.write y 4 ((0.5)*((V.!) x' 4))
  MV.write y 5 ((0.25)*((V.!) x' 5))
  MV.write y 6 ((0.125)*((V.!) x' 6))
  
  return ()

diagMatLMCheck :: Bool
diagMatLMCheck = let (converged, ePairs) = unsafePerformIO $eigs diagMat 1000 LM 6 (0) 6000
                     evals = map fst ePairs
                     evecs = map snd ePairs
                     trueEvals = map (:+0) [1000,999..]
                     diff = map (realPart . abs) $ zipWith (-) evals trueEvals
                 in not $ and (map (>=1e-10) diff)

diagMatSMCheck :: Bool
diagMatSMCheck = let (converged, ePairs) = unsafePerformIO $eigs diagMat 1000 SM 6 (0) 6000
                     evals = map fst ePairs
                     evecs = map snd ePairs
                     trueEvals = map (:+0) [1,2..]
                     diff = map (realPart . abs) $ zipWith (-) evals trueEvals
                 in not $ and (map (>=1e-10) diff)

--complexEigsCheck :: Bool
complexEigsCheck = let (converged, ePairs) = unsafePerformIO $ eigs complexEigs 7 LM 5 (0) 70
                       evals = map fst ePairs
                       evecs = map snd ePairs
                       trueEvals = [1 :+ 2, 1 :+ (-2), 2, 1, 0.5]
                       diff = map (realPart . abs) $ zipWith (-) evals trueEvals
                   in not $ and (map (>=1e-10) diff) 

--not $ and (map (>=1e-10) diff)