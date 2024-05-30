# CSE305 - Project

To run the code you must run the following commands:

`tar -xzvf ImageMagick.tar.gz`

`cd ImageMagick-7.*`

`./configure --prefix=$HOME/ImageMagick`

`make`

`make install`

`export PATH=$HOME/ImageMagick/bin:$PATH`

`export PKG_CONFIG_PATH=$HOME/ImageMagick/lib/pkgconfig:$PKG_CONFIG_PATH`

`export LD_LIBRARY_PATH=$HOME/ImageMagick/lib:$LD_LIBRARY_PATH`

and then run

`find $HOME/ImageMagick -name "Magick++.h"`

copy the result in `<found_include_path>`:

`g++ -std=c++11 -o main main.cpp -I<found_include_path> -L$HOME/ImageMagick/lib -lMagick++-7.Q16HDRI -lMagickCore-7.Q16HDRI -lMagickWand-7.Q16HDRI`

and then

`./main`