{-# LANGUAGE ForeignFunctionInterface #-}

module Numeric.LinearAlgebra.UnsymmetricEigen (eigs, 
                                               Which(LM, SM),  
                                               ArpackLinearOp) where

import Foreign
import Foreign.Ptr
import Foreign.Marshal.Array
import Foreign.C.Types
import Foreign.C.String
import Foreign.Storable
import Foreign.Storable.Complex
import Foreign.ForeignPtr
import qualified Data.Vector.Storable as V
import qualified Data.Vector.Storable.Mutable as SV
import Control.Monad.Loops
import Control.Monad
import Data.Complex
 
data Which = LM | SM deriving Show

data ArpackSetup = ArpackSetup { ido :: (Ptr CInt),
                                 bmat :: (Ptr CChar),
                                 n :: (Ptr CInt),
                                 which :: (Ptr CChar),
                                 nev :: (Ptr CInt),
                                 tol :: (Ptr CDouble),
                                 resid :: (Ptr CDouble),
                                 ncv :: (Ptr CInt),
                                 v :: (Ptr CDouble),
                                 ldv :: (Ptr CInt),
                                 iparam :: (Ptr CInt),
                                 ipntr :: (Ptr CInt),
                                 workd :: (Ptr CDouble),
                                 workl :: (Ptr CDouble),
                                 lworkl :: (Ptr CInt),
                                 info :: (Ptr CInt) } deriving (Show)
                                                               
type ProblemDim = Int
type NumEV = Int
type Tolerance = Double
type MaxIter = Int

setUpArpack :: ProblemDim -> Which -> NumEV -> Tolerance -> MaxIter -> 
               IO ArpackSetup
setUpArpack n which nev tol maxIter = do
  let ncv = min n (2*nev + 1)
      lWorkl = 3*ncv^2 + 6*ncv
      ldv = n
  [idoPtr,nPtr,nevPtr,ncvPtr,ldvPtr] <- sequence $ replicate 5 malloc :: IO [Ptr CInt]
  [bmatPtr, whichPtr] <- sequence [mallocArray 2, mallocArray 3] :: IO [Ptr CChar]
  [tolPtr, residPtr, vPtr, workdPtr, worklPtr] <- 
    sequence $ map mallocArray [1, n, (n*ncv), (3*n), lWorkl] :: IO [Ptr CDouble]
  [iParamPtr, iPntrPtr, lworklPtr, infoPtr] <- 
    sequence $ map mallocArray [11, 14, 1, 1] :: IO [Ptr CInt]
  
  poke idoPtr 0
  poke bmatPtr $ castCharToCChar 'I'
  poke nPtr $ fromIntegral n
  pokeArray whichPtr (map castCharToCChar (show which))
  poke nevPtr $ fromIntegral nev
  poke tolPtr $ realToFrac tol
  poke ncvPtr $ fromIntegral ncv
  poke ldvPtr $ fromIntegral ldv
  poke iParamPtr $ fromIntegral 1
  pokeElemOff iParamPtr 2 $ fromIntegral maxIter
  pokeElemOff iParamPtr 3 $ fromIntegral 1
  pokeElemOff iParamPtr 6 $ fromIntegral 1
  poke lworklPtr $ fromIntegral lWorkl
  poke infoPtr $ fromIntegral 0
  
  c_dnaupd idoPtr bmatPtr nPtr whichPtr nevPtr tolPtr residPtr ncvPtr vPtr
    ldvPtr iParamPtr iPntrPtr workdPtr worklPtr lworklPtr infoPtr
  
  return (ArpackSetup idoPtr bmatPtr nPtr whichPtr nevPtr tolPtr residPtr ncvPtr vPtr
    ldvPtr iParamPtr iPntrPtr workdPtr worklPtr lworklPtr infoPtr)

type ArpackLinearOp = (SV.IOVector CDouble -> SV.IOVector CDouble -> IO ())


iterateArpack :: ArpackLinearOp -> ArpackSetup -> IO ()  
iterateArpack f ar = do
  workdXElem <- peek (ipntr ar)
  workdYElem <- peekElemOff (ipntr ar) 1
  workdZElem <- peekElemOff (ipntr ar) 2
  n' <- peek (n ar)
  xPtr <- newForeignPtr_ $ advancePtr (workd ar) ((fromIntegral workdXElem) - 1)
  yPtr <- newForeignPtr_ $ advancePtr (workd ar) ((fromIntegral workdYElem) - 1)
  --zPtr <- newForeignPtr_ $ advancePtr (workd ar) ((fromIntegral workdZElem) - 1)
  let x = SV.unsafeFromForeignPtr xPtr 0 (fromIntegral n')
      y = SV.unsafeFromForeignPtr yPtr 0 (fromIntegral n')
      --z = SV.unsafeFromForeignPtr zPtr 0 (fromIntegral n')
      idoPtr = ido ar 
      bmatPtr = bmat ar
      nPtr = n ar
      whichPtr = which ar
      nevPtr = nev ar 
      tolPtr = tol ar
      residPtr = resid ar
      ncvPtr = ncv ar
      vPtr = v ar
      ldvPtr = ldv ar
      iParamPtr = iparam ar
      iPntrPtr = ipntr ar
      workdPtr = workd ar
      worklPtr = workl ar
      lworklPtr = lworkl ar
      infoPtr = info ar
  --SV.copy z x
  f x y
   
  c_dnaupd idoPtr bmatPtr nPtr whichPtr nevPtr tolPtr residPtr ncvPtr vPtr
    ldvPtr iParamPtr iPntrPtr workdPtr worklPtr lworklPtr infoPtr
  
  return ()

arpack :: ArpackLinearOp -> ProblemDim -> Which -> NumEV -> Tolerance -> 
          MaxIter -> IO ArpackSetup
arpack f n which nev tol mxItr = do
  ar <- setUpArpack n which nev tol mxItr
  whileM_ ((peek $ ido ar) >>= (\ido' -> return (ido' == 1))) 
    $ iterateArpack f ar
  return ar
  
errCheck :: ArpackResults -> IO Bool
errCheck ar = do 
  info' <- (peek $ info $ arSetup ar)
  putStrLn $ show info'
  return (info' < 0)

data ArpackResults = ArpackResults { rvec :: Ptr CInt,
                                     howmny :: Ptr CChar,
                                     select :: Ptr CInt,
                                     dr :: SV.IOVector CDouble,
                                     di :: SV.IOVector CDouble, 
                                     z :: SV.IOVector CDouble,
                                     ldz :: Ptr CInt,
                                     sigmar :: Ptr CDouble,
                                     sigmai :: Ptr CDouble,
                                     workev :: Ptr CDouble,
                                     arSetup :: ArpackSetup}


parseArpackOutput :: ArpackSetup -> IO ArpackResults 
parseArpackOutput ar = do
  ncv' <- peek $ ncv ar
  nev' <- peek $ nev ar
  n' <- peek $ n ar
  
  rvecPtr <- malloc :: IO (Ptr CInt)
  howmnyPtr <- malloc :: IO (Ptr CChar)
  drVec <- SV.new ((fromIntegral nev')) :: IO (SV.IOVector CDouble)   
  diVec <- SV.new ((fromIntegral nev')) :: IO (SV.IOVector CDouble)
  selectPtr <- mallocArray (fromIntegral ncv') :: IO (Ptr CInt)
  zVec <- SV.new (fromIntegral (n'*(nev'+1))) :: IO (SV.IOVector CDouble)
  ldzPtr <- malloc :: IO (Ptr CInt)
  sigmarPtr <- malloc :: IO (Ptr CDouble)
  sigmaiPtr <- malloc :: IO (Ptr CDouble)
  workevPtr <- mallocArray (fromIntegral (3*ncv')) :: IO (Ptr CDouble)
  
  poke rvecPtr 1
  poke howmnyPtr $ castCharToCChar 'A'
  poke ldzPtr n'
  
  let (drfPtr,_,_) = SV.unsafeToForeignPtr drVec
      drPtr = unsafeForeignPtrToPtr drfPtr
      (difPtr,_,_) = SV.unsafeToForeignPtr diVec
      diPtr = unsafeForeignPtrToPtr difPtr
      (zfPtr,_,_) = SV.unsafeToForeignPtr zVec
      zPtr = unsafeForeignPtrToPtr zfPtr
  c_dneupd rvecPtr howmnyPtr selectPtr drPtr diPtr zPtr ldzPtr sigmarPtr
    sigmaiPtr workevPtr (bmat ar) (n ar) (which ar) (nev ar) (tol ar)
    (resid ar) (ncv ar) (v ar) (ldv ar) (iparam ar) (ipntr ar) (workd ar) 
    (workl ar) (lworkl ar) (info ar)
  
  return $ ArpackResults rvecPtr howmnyPtr selectPtr drVec diVec zVec ldzPtr
    sigmarPtr sigmaiPtr workevPtr ar
  
getEigenValue :: ArpackResults -> Int -> IO (Complex Double)
getEigenValue ar idx = do
  let reEigs = dr ar
      imEigs = di ar
  re <- SV.read reEigs idx
  im <- SV.read imEigs idx
  return ((realToFrac re :: Double) :+ (realToFrac im :: Double))


getEigenVectorReal :: ArpackResults -> Int -> IO (V.Vector (Complex Double))
getEigenVectorReal ar idx = do
  n' <- peek $ n $ arSetup ar
  ldz' <- peek $ ldz ar
  let z' = z ar
  evec <- V.freeze $ SV.slice (idx*(fromIntegral ldz')) (fromIntegral n') z'
  return $ V.map ((:+0) . realToFrac) evec



getEigenVectorComplex :: ArpackResults -> Int -> Bool -> IO (V.Vector (Complex Double))
getEigenVectorComplex ar idx conj = do
  n' <- peek $ n $ arSetup ar
  ldz' <- peek $ ldz ar
  let z' = z ar
  reEvec <- V.freeze $ SV.slice (idx*(fromIntegral ldz')) (fromIntegral n') z'
  imEvec <- V.freeze $ SV.slice ((idx+1)*(fromIntegral ldz')) (fromIntegral n') z'
  if conj
     then return $ V.zipWith (:+) (V.map realToFrac reEvec) (V.map ((*(-1)) . realToFrac) imEvec)
    else return $ V.zipWith (:+) (V.map realToFrac reEvec) (V.map realToFrac imEvec)
  

getEigenPair :: ArpackResults -> Int -> IO (Complex Double, V.Vector (Complex Double))
getEigenPair ar idx = do 
  eig <- getEigenValue ar idx
  let evec = case (compare (imagPart eig) 0) of 
        EQ -> getEigenVectorReal ar idx 
        LT -> getEigenVectorComplex ar (idx-1) True
        GT -> getEigenVectorComplex ar idx False
  evec' <- evec
  return (eig,evec')

eigs :: ArpackLinearOp -> ProblemDim -> Which -> NumEV -> Tolerance -> 
        MaxIter -> IO (Bool, [(Complex Double, V.Vector (Complex Double))])
eigs f n which nev tol' iters = do
  ar <- arpack f n which nev tol' iters
  arOut <- parseArpackOutput ar
  err <- errCheck arOut
  if (err)
    then do return (False,[]) 
    else do
      pairs <- sequence $ map (getEigenPair arOut) [0..(nev-1)] 
      return (True,pairs)

foreign import ccall unsafe "arpack.h dnaupd_"
  c_dnaupd :: Ptr CInt -> Ptr CChar -> Ptr CInt -> Ptr CChar -> Ptr CInt ->
              Ptr CDouble -> Ptr CDouble -> Ptr CInt -> Ptr CDouble ->
              Ptr CInt -> Ptr CInt -> Ptr CInt -> Ptr CDouble ->
              Ptr CDouble -> Ptr CInt -> Ptr CInt -> IO ()
                 
foreign import ccall unsafe "arpack.h dneupd_"
  c_dneupd :: Ptr CInt -> Ptr CChar -> Ptr CInt -> Ptr CDouble -> Ptr CDouble -> 
              Ptr CDouble -> Ptr CInt -> Ptr CDouble -> Ptr CDouble -> 
              Ptr CDouble -> Ptr CChar -> Ptr CInt -> Ptr CChar -> Ptr CInt -> 
              Ptr CDouble -> Ptr CDouble -> Ptr CInt -> Ptr CDouble -> 
              Ptr CInt -> Ptr CInt -> Ptr CInt -> Ptr CDouble -> Ptr CDouble -> 
              Ptr CInt -> Ptr CInt -> IO ()
              
