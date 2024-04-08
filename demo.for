C This is a simple multiplication program in Fortran 77
      PROGRAM MULTIPLY
      INTEGER A, B, RESULT
      PRINT *, 'Enter two numbers to multiply:'
      READ (*,*) A, B
      RESULT = A * B
      PRINT *, 'The result of the multiplication is:', RESULT
      END