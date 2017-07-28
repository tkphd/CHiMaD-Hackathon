#!/bin/bash

if [[ -f ../square/electrochem.000100.dat ]]
then
	mmsp2png --field=0 ../square/electrochem.000100.dat square_c.005.png
	mmsp2png --field=2 ../square/electrochem.000100.dat square_p.005.png
fi

if [[ -f ../square/electrochem.000200.dat ]]
then
	mmsp2png --field=0 ../square/electrochem.000200.dat square_c.010.png
	mmsp2png --field=2 ../square/electrochem.000200.dat square_p.010.png
fi

if [[ -f ../square/electrochem.000400.dat ]]
then
	mmsp2png --field=0 ../square/electrochem.000400.dat square_c.020.png
	mmsp2png --field=2 ../square/electrochem.000400.dat square_p.020.png
fi

if [[ -f ../square/electrochem.001000.dat ]]
then
	mmsp2png --field=0 ../square/electrochem.001000.dat square_c.050.png
	mmsp2png --field=2 ../square/electrochem.001000.dat square_p.050.png
fi

if [[ -f ../square/electrochem.002000.dat ]]
then
	mmsp2png --field=0 ../square/electrochem.002000.dat square_c.100.png
	mmsp2png --field=2 ../square/electrochem.002000.dat square_p.100.png
fi

if [[ -f ../square/electrochem.004000.dat ]]
then
	mmsp2png --field=0 ../square/electrochem.004000.dat square_c.200.png
	mmsp2png --field=2 ../square/electrochem.004000.dat square_p.200.png
fi

if [[ -f ../square/electrochem.008000.dat ]]
then
	mmsp2png --field=0 ../square/electrochem.008000.dat square_c.400.png
	mmsp2png --field=2 ../square/electrochem.008000.dat square_p.400.png
fi

if [[ -f ../square/electrochem.020000.dat ]]
then
	mmsp2png --field=0 ../square/electrochem.020000.dat square_c.1000.png
	mmsp2png --field=2 ../square/electrochem.020000.dat square_p.1000.png
fi
