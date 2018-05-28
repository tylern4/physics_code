*
* $Id: gammf.F,v 1.1.1.1 1996/04/01 15:01:54 mclareni Exp $
*
* $Log: gammf.F,v $
* Revision 1.1.1.1  1996/04/01 15:01:54  mclareni
* Mathlib gen
*
*
      FUNCTION GAMMF(X)
      DOUBLE PRECISION DGAMMF,D
      SROUND(D)=D+(D-SNGL(D))
      GAMMF=SROUND(DGAMMF(DBLE(X)))
      RETURN
      END
