
#include <stdio.h>
#include <iostream>
#include <vector>

#include <colormapper.h>

static irwincolor::Map c;

extern "C" int load_colormap_(char* cfile, char* cmapname)
{
	int io = 0;

	std::string file = cfile;
	std::string mapname = cmapname;
	if (file != "" && mapname != "") io = c.load(file, mapname);

	return io;
}

extern "C" void map_(double& x, uint8_t* rgbo)
{
	//std::cout << "x = " << x << std::endl;

	std::vector<uint8_t> rgb = c.map(x);

	rgbo[0] = rgb[0];
	rgbo[1] = rgb[1];
	rgbo[2] = rgb[2];

	return;
}

extern "C" int writepng_(uint8_t* b, int& nx, int& ny, char* cf)
{
	// Can pix be constructed without copying every element of b?
	// May need to add alpha channel in Fortran.

	//std::vector<uint8_t> pix(b, b + 4 * nx * ny);
	std::vector<uint8_t> pix(4 * nx * ny);
	int j = 0;
	for (int i = 0; i < 4 * nx * ny; i++)
	{
		if ((i+1) % 4 != 0)
			pix[i] = b[j++];
			// or
			//   pix[j * 4/3] = b[j];
		else
			pix[i] = 255;
	}

	return irwincolor::savePng(pix, nx, ny, (std::string) cf);
}

extern "C" void sortidx_(int* v, int& n, int* idx)
{
	// Wrap the C++ sortidx function, which takes/returns vectors, with this
	// function that takes/returns simple Fortran-compatible arrays.
	//
	// There is overhead for copying array/vectors, so it would be better to
	// implement sorting natively in Fortran.

	//std::vector<int> vl(std::begin(*v), std::end(*v));
	std::vector<int> vl(v, v + n);

	//std::cout << "vl = ";
	//for (int i = 0; i < 100; i++)
	//	std::cout << vl[i] << " ";
	//std::cout << std::endl;

	auto idxl = irwincolor::sortidx(vl);

	//std::cout << "idxl = ";
	//for (int i = 0; i < 100; i++)
	//	std::cout << idxl[i] << " ";
	//std::cout << std::endl;

	std::copy(idxl.begin(), idxl.end(), idx);

}

