reset
set grid 

pl    'Re10000/taxafft.dat' u 1:2 w lp,  \
      'Re15000/taxafft.dat' u 1:2 w lp,  \
      'Re20000/taxafft.dat' u 1:2 w lp,  \
      '../case_test_euler/taxafft.dat' u 1:2 w lp

pause mouse keypress "Type a letter from A-F in the active window"
