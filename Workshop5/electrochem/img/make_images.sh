#!/bin/bash

# Process files with dx = 1.0

if [[ -f ../square/electrochem_dx10.00100.dat ]]
then
	if [[ ! -f square_dx10_c.0005.png ]]
	then
		mmsp2png --field=0 ../square/electrochem_dx10.00100.dat square_dx10_c.0005.png
		mmsp2png --field=2 ../square/electrochem_dx10.00100.dat square_dx10_p.0005.png
	fi
	if [[ ! -f ../square/electrochem_dx10.00100.csv ]]
	then
		.././linescan ../square/electrochem_dx10.00100.dat ../square/electrochem_dx10.00100.csv
	fi
fi

if [[ -f ../square/electrochem_dx10.00200.dat ]]
then
	if [[ ! -f square_dx10_c.0010.png ]]
	then
		mmsp2png --field=0 ../square/electrochem_dx10.00200.dat square_dx10_c.0010.png
		mmsp2png --field=2 ../square/electrochem_dx10.00200.dat square_dx10_p.0010.png
	fi
	if [[ ! -f ../square/electrochem_dx10.00200.csv ]]
	then
		.././linescan ../square/electrochem_dx10.00200.dat ../square/electrochem_dx10.00200.csv
	fi
fi

if [[ -f ../square/electrochem_dx10.00400.dat ]]
then
	if [[ ! -f square_dx10_c.0020.png ]]
	then
		mmsp2png --field=0 ../square/electrochem_dx10.00400.dat square_dx10_c.0020.png
		mmsp2png --field=2 ../square/electrochem_dx10.00400.dat square_dx10_p.0020.png
	fi
	if [[ ! -f ../square/electrochem_dx10.00400.csv ]]
	then
		.././linescan ../square/electrochem_dx10.00400.dat ../square/electrochem_dx10.00400.csv
	fi
fi

if [[ -f ../square/electrochem_dx10.01000.dat ]]
then
	if [[ ! -f square_dx10_c.0050.png ]]
	then
		mmsp2png --field=0 ../square/electrochem_dx10.01000.dat square_dx10_c.0050.png
		mmsp2png --field=2 ../square/electrochem_dx10.01000.dat square_dx10_p.0050.png
	fi
	if [[ ! -f ../square/electrochem_dx10.01000.csv ]]
	then
		.././linescan ../square/electrochem_dx10.01000.dat ../square/electrochem_dx10.01000.csv
	fi
fi

if [[ -f ../square/electrochem_dx10.02000.dat ]]
then
	if [[ ! -f square_dx10_c.0100.png ]]
	then
		mmsp2png --field=0 ../square/electrochem_dx10.02000.dat square_dx10_c.0100.png
		mmsp2png --field=2 ../square/electrochem_dx10.02000.dat square_dx10_p.0100.png
	fi
	if [[ ! -f ../square/electrochem_dx10.02000.csv ]]
	then
		.././linescan ../square/electrochem_dx10.02000.dat ../square/electrochem_dx10.02000.csv
	fi
fi

if [[ -f ../square/electrochem_dx10.04000.dat ]]
then
	if [[ ! -f square_dx10_c.0200.png ]]
	then
		mmsp2png --field=0 ../square/electrochem_dx10.04000.dat square_dx10_c.0200.png
		mmsp2png --field=2 ../square/electrochem_dx10.04000.dat square_dx10_p.0200.png
	fi
	if [[ ! -f ../square/electrochem_dx10.04000.csv ]]
	then
		.././linescan ../square/electrochem_dx10.04000.dat ../square/electrochem_dx10.04000.csv
	fi
fi

if [[ -f ../square/electrochem_dx10.08000.dat ]]
then
	if [[ ! -f square_dx10_c.0400.png ]]
	then
		mmsp2png --field=0 ../square/electrochem_dx10.08000.dat square_dx10_c.0400.png
		mmsp2png --field=2 ../square/electrochem_dx10.08000.dat square_dx10_p.0400.png
	fi
	if [[ ! -f ../square/electrochem_dx10.08000.csv ]]
	then
		.././linescan ../square/electrochem_dx10.08000.dat ../square/electrochem_dx10.08000.csv
	fi
fi

if [[ -f ../square/electrochem_dx10.20000.dat ]]
then
	if [[ ! -f square_dx10_c.1000.png ]]
	then
		mmsp2png --field=0 ../square/electrochem_dx10.20000.dat square_dx10_c.1000.png
		mmsp2png --field=2 ../square/electrochem_dx10.20000.dat square_dx10_p.1000.png
	fi
	if [[ ! -f ../square/electrochem_dx10.20000.csv ]]
	then
		.././linescan ../square/electrochem_dx10.20000.dat ../square/electrochem_dx10.20000.csv
	fi
fi

# Process files with dx = 0.5

if [[ -f ../square/electrochem_dx05.00100.dat ]]
then
	if [[ ! -f square_dx05_c.0005.png ]]
	then
		mmsp2png --field=0 ../square/electrochem_dx05.00100.dat square_dx05_c.0005.png
		mmsp2png --field=2 ../square/electrochem_dx05.00100.dat square_dx05_p.0005.png
	fi
	if [[ ! -f ../square/electrochem_dx05.00100.csv ]]
	then
		.././linescan ../square/electrochem_dx05.00100.dat ../square/electrochem_dx05.00100.csv
	fi
fi

if [[ -f ../square/electrochem_dx05.00200.dat ]]
then
	if [[ ! -f square_dx05_c.0010.png ]]
	then
		mmsp2png --field=0 ../square/electrochem_dx05.00200.dat square_dx05_c.0010.png
		mmsp2png --field=2 ../square/electrochem_dx05.00200.dat square_dx05_p.0010.png
	fi
	if [[ ! -f ../square/electrochem_dx05.00200.csv ]]
	then
		.././linescan ../square/electrochem_dx05.00200.dat ../square/electrochem_dx05.00200.csv
	fi
fi

if [[ -f ../square/electrochem_dx05.00400.dat ]]
then
	if [[ ! -f square_dx05_c.0020.png ]]
	then
		mmsp2png --field=0 ../square/electrochem_dx05.00400.dat square_dx05_c.0020.png
		mmsp2png --field=2 ../square/electrochem_dx05.00400.dat square_dx05_p.0020.png
	fi
	if [[ ! -f ../square/electrochem_dx05.00400.csv ]]
	then
		.././linescan ../square/electrochem_dx05.00400.dat ../square/electrochem_dx05.00400.csv
	fi
fi

if [[ -f ../square/electrochem_dx05.01000.dat ]]
then
	if [[ ! -f square_dx05_c.0050.png ]]
	then
		mmsp2png --field=0 ../square/electrochem_dx05.01000.dat square_dx05_c.0050.png
		mmsp2png --field=2 ../square/electrochem_dx05.01000.dat square_dx05_p.0050.png
	fi
	if [[ ! -f ../square/electrochem_dx05.01000.csv ]]
	then
		.././linescan ../square/electrochem_dx05.01000.dat ../square/electrochem_dx05.01000.csv
	fi
fi

if [[ -f ../square/electrochem_dx05.02000.dat ]]
then
	if [[ ! -f square_dx05_c.0100.png ]]
	then
		mmsp2png --field=0 ../square/electrochem_dx05.02000.dat square_dx05_c.0100.png
		mmsp2png --field=2 ../square/electrochem_dx05.02000.dat square_dx05_p.0100.png
	fi
	if [[ ! -f ../square/electrochem_dx05.02000.csv ]]
	then
		.././linescan ../square/electrochem_dx05.02000.dat ../square/electrochem_dx05.02000.csv
	fi
fi

if [[ -f ../square/electrochem_dx05.04000.dat ]]
then
	if [[ ! -f square_dx05_c.0200.png ]]
	then
		mmsp2png --field=0 ../square/electrochem_dx05.04000.dat square_dx05_c.0200.png
		mmsp2png --field=2 ../square/electrochem_dx05.04000.dat square_dx05_p.0200.png
	fi
	if [[ ! -f ../square/electrochem_dx05.04000.csv ]]
	then
		.././linescan ../square/electrochem_dx05.04000.dat ../square/electrochem_dx05.04000.csv
	fi
fi

if [[ -f ../square/electrochem_dx05.08000.dat ]]
then
	if [[ ! -f square_dx05_c.0400.png ]]
	then
		mmsp2png --field=0 ../square/electrochem_dx05.08000.dat square_dx05_c.0400.png
		mmsp2png --field=2 ../square/electrochem_dx05.08000.dat square_dx05_p.0400.png
	fi
	if [[ ! -f ../square/electrochem_dx05.08000.csv ]]
	then
		.././linescan ../square/electrochem_dx05.08000.dat ../square/electrochem_dx05.08000.csv
	fi
fi

if [[ -f ../square/electrochem_dx05.20000.dat ]]
then
	if [[ ! -f square_dx05_c.1000.png ]]
	then
		mmsp2png --field=0 ../square/electrochem_dx05.20000.dat square_dx05_c.1000.png
		mmsp2png --field=2 ../square/electrochem_dx05.20000.dat square_dx05_p.1000.png
	fi
	if [[ ! -f ../square/electrochem_dx05.20000.csv ]]
	then
		.././linescan ../square/electrochem_dx05.20000.dat ../square/electrochem_dx05.20000.csv
	fi
fi

# Process files with dx = 0.1

if [[ -f ../square/electrochem_dx01.00100.dat ]]
then
	if [[ ! -f square_dx01_c.0005.png ]]
	then
		mmsp2png --field=0 ../square/electrochem_dx01.00100.dat square_dx01_c.0005.png
		mmsp2png --field=2 ../square/electrochem_dx01.00100.dat square_dx01_p.0005.png
	fi
	if [[ ! -f ../square/electrochem_dx01.00100.csv ]]
	then
		.././linescan ../square/electrochem_dx01.00100.dat ../square/electrochem_dx01.00100.csv
	fi
fi

if [[ -f ../square/electrochem_dx01.00200.dat ]]
then
	if [[ ! -f square_dx01_c.0010.png ]]
	then
		mmsp2png --field=0 ../square/electrochem_dx01.00200.dat square_dx01_c.0010.png
		mmsp2png --field=2 ../square/electrochem_dx01.00200.dat square_dx01_p.0010.png
	fi
	if [[ ! -f ../square/electrochem_dx01.00200.csv ]]
	then
		.././linescan ../square/electrochem_dx01.00200.dat ../square/electrochem_dx01.00200.csv
	fi
fi

if [[ -f ../square/electrochem_dx01.00400.dat ]]
then
	if [[ ! -f square_dx01_c.0020.png ]]
	then
		mmsp2png --field=0 ../square/electrochem_dx01.00400.dat square_dx01_c.0020.png
		mmsp2png --field=2 ../square/electrochem_dx01.00400.dat square_dx01_p.0020.png
	fi
	if [[ ! -f ../square/electrochem_dx01.00400.csv ]]
	then
		.././linescan ../square/electrochem_dx01.00400.dat ../square/electrochem_dx01.00400.csv
	fi
fi

if [[ -f ../square/electrochem_dx01.01000.dat ]]
then
	if [[ ! -f square_dx01_c.0050.png ]]
	then
		mmsp2png --field=0 ../square/electrochem_dx01.01000.dat square_dx01_c.0050.png
		mmsp2png --field=2 ../square/electrochem_dx01.01000.dat square_dx01_p.0050.png
	fi
	if [[ ! -f ../square/electrochem_dx01.01000.csv ]]
	then
		.././linescan ../square/electrochem_dx01.01000.dat ../square/electrochem_dx01.01000.csv
	fi
fi

if [[ -f ../square/electrochem_dx01.02000.dat ]]
then
	if [[ ! -f square_dx01_c.0100.png ]]
	then
		mmsp2png --field=0 ../square/electrochem_dx01.02000.dat square_dx01_c.0100.png
		mmsp2png --field=2 ../square/electrochem_dx01.02000.dat square_dx01_p.0100.png
	fi
	if [[ ! -f ../square/electrochem_dx01.02000.csv ]]
	then
		.././linescan ../square/electrochem_dx01.02000.dat ../square/electrochem_dx01.02000.csv
	fi
fi

if [[ -f ../square/electrochem_dx01.04000.dat ]]
then
	if [[ ! -f square_dx01_c.0200.png ]]
	then
		mmsp2png --field=0 ../square/electrochem_dx01.04000.dat square_dx01_c.0200.png
		mmsp2png --field=2 ../square/electrochem_dx01.04000.dat square_dx01_p.0200.png
	fi
	if [[ ! -f ../square/electrochem_dx01.04000.csv ]]
	then
		.././linescan ../square/electrochem_dx01.04000.dat ../square/electrochem_dx01.04000.csv
	fi
fi

if [[ -f ../square/electrochem_dx01.08000.dat ]]
then
	if [[ ! -f square_dx01_c.0400.png ]]
	then
		mmsp2png --field=0 ../square/electrochem_dx01.08000.dat square_dx01_c.0400.png
		mmsp2png --field=2 ../square/electrochem_dx01.08000.dat square_dx01_p.0400.png
	fi
	if [[ ! -f ../square/electrochem_dx01.08000.csv ]]
	then
		.././linescan ../square/electrochem_dx01.08000.dat ../square/electrochem_dx01.08000.csv
	fi
fi

if [[ -f ../square/electrochem_dx01.20000.dat ]]
then
	if [[ ! -f square_dx01_c.1000.png ]]
	then
		mmsp2png --field=0 ../square/electrochem_dx01.20000.dat square_dx01_c.1000.png
		mmsp2png --field=2 ../square/electrochem_dx01.20000.dat square_dx01_p.1000.png
	fi
	if [[ ! -f ../square/electrochem_dx01.20000.csv ]]
	then
		.././linescan ../square/electrochem_dx01.20000.dat ../square/electrochem_dx01.20000.csv
	fi
fi
