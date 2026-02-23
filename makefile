# Copyright 2025 Giorgio Gamba

# NOTE: $@ stands for the files to build (main.cpp), $^ is for its dependencies

# Configuration variables to be edited as you prefer
BuildPath =  build
ExecutableName = avil

LibraryFlags = -I ./lib/portaudio/include\
				  ./lib/portaudio/lib/.libs/libportaudio.a -framework CoreAudio\
				  -framework AudioToolbox -framework CoreServices -pthread\
				  -I ./lib/fftw-3.3.10/api -lfftw3\
				  -I ./lib/libsndfile-1.2.2/include -lm -L ./lib/libsndfile-1.2.2/src/.libs -lsndfile\
				  -I ./lib/FTXUI-5.0.0/include -L ./lib/FTXUI-5.0.0/build/build/lib -lftxui-component -lftxui-dom -lftxui-screen

$(BuildPath)/$(ExecutableName): src/main.cpp src/fft.hpp src/types.h src/constants.h
	g++ -o $@ src/main.cpp $(LibraryFlags) --std=c++20 -mavx2 -mfma
	
# By setting this command, all the targets will be executed consequently
install_dependencies: install_portaudio install_fftw install_libsndfile install_ftxui
.PHONY: install_dependencies

uninstall_dependencies: uninstall_portaudio uninstall_fftw uninstall_libsndfile uninstall_ftxui
.PHONY: uninstall_dependencies

install_portaudio:
	mkdir -p lib 
	curl https://files.portaudio.com/archives/pa_stable_v190700_20210406.tgz | tar -zx -C lib
	cd lib/portaudio && ./configure && $(MAKE) -j
.PHONY: install_portaudio

uninstall_portaudio:
	cd lib/portaudio && $(MAKE) uninstall
	rm -rf lib/portaudio
.PHONY: uninstall_portaudio

install_fftw:
	mkdir -p lib
	curl https://www.fftw.org/fftw-3.3.10.tar.gz | tar -zx -C lib
	cd lib/fftw-3.3.10 && ./configure && $(MAKE) -j && sudo $(MAKE) install
.PHONY: install_fftw 

uninstall_fftw:
	cd lib/fftw-3.3.10 && $(MAKE) uninstall
	rm -rf lib/fftw-3.3.10
.PHONY: uninstall_fftw

install_libsndfile:
	mkdir -p lib
	curl -L https://github.com/libsndfile/libsndfile/releases/download/1.2.2/libsndfile-1.2.2.tar.xz | tar -Jx -C lib
	cd lib/libsndfile-1.2.2 && ./configure && $(MAKE)
.PHONY: install_libsndfile

uninstall_libsndfile:
	cd lib/libsndfile-1.2.2 && $(MAKE) uninstall
	rm -rf lib/libsndfile-1.2.2
.PHONY: uninstall_libsndfile

install_ftxui:
	mkdir -p lib
	curl -L https://github.com/ArthurSonzogni/FTXUI/archive/refs/tags/v5.0.0.tar.gz | tar -zx -C lib
	cd lib/FTXUI-5.0.0 && cmake -S . -B build \
	-DCMAKE_BUILD_TYPE=Release \
	-DCMAKE_ARCHIVE_OUTPUT_DIRECTORY=build/lib \
	-DFTXUI_BUILD_EXAMPLES=OFF \
	-DFTXUI_BUILD_DOCS=OFF \
	-DFTXUI_BUILD_TESTS=OFF && \
	cmake --build build
.PHONY: install_ftxui

uninstall_ftxui:
	rm -rf lib/FTXUI-5.0.0
.PHONY: uninstall_libsndfile

clean:
	rm -f $(BuildPath)/$(ExecutableName)
.PHONY: clean
