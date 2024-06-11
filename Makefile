PREFIX=$(HOME)/ImageMagick
IM_VERSION=7.*

all: configure build install setenv

configure:
	cd ImageMagick-$(IM_VERSION) && ./configure --prefix=$(PREFIX)

build:
	cd ImageMagick-$(IM_VERSION) && make

install:
	cd ImageMagick-$(IM_VERSION) && make install

setenv:
	@echo "export PATH=$(PREFIX)/bin:$$PATH" >> ~/.bashrc
	@echo "export PKG_CONFIG_PATH=$(PREFIX)/lib/pkgconfig:$$PKG_CONFIG_PATH" >> ~/.bashrc
	@echo "export LD_LIBRARY_PATH=$(PREFIX)/lib:$$LD_LIBRARY_PATH" >> ~/.bashrc

clean:
	cd ImageMagick-$(IM_VERSION) && make clean

.PHONY: all configure build install setenv clean