#!/bin/bash
echo "Running tests..."

# Esegui i test
echo " testing the fortran library..."
./test_fortran/test_fortran.exe 
diff -q output_fortran.txt ref/output_fortran.ref > output_fortran.diff
F_RESULT=$?
if [ $F_RESULT -eq 0 ]; then
  echo "✅ fortran test passed."
fi
echo
echo " testing the c interface..."
./test_c/test_c.exe 

# Confronta con i file di riferimento
diff -q output_c.txt ref/output_c.ref > output_c.diff
C_RESULT=$?
if [ $C_RESULT -eq 0 ]; then
  echo "✅ C interface test passed."
fi
echo

if [ $C_RESULT -eq 0 ] && [ $F_RESULT -eq 0 ]; then
  echo "✅ All tests are passing."
  rm output_fortran.* output_c.*
  exit 0
else
  echo "❌ Failed Tests:"
  [ $C_RESULT -ne 0 ] && echo "  - test_c "
  [ $F_RESULT -ne 0 ] && echo "  - test_fortran "
  exit 1
fi
