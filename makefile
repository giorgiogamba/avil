# Copyright 2025 Giorgio Gamba

# NOTE: $@ stands for the files to build (main.cpp), $^ is for its dependencies

# Configuration variables to be edited as you prefer
BuildPath =  build
ExecutableName = avil

LibraryFlags = -I ./lib/portaudio/include\
				  ./lib/portaudio/lib/.libs/libportaudio.a -framework CoreAudio\
				  -framework AudioToolbox -framework CoreServices -pthread\
				  -I ./lib/fftw-3.3.10/api -lfftw3

$(BuildPath)/$(ExecutableName): src/main.cpp
	g++ -o $@ $^ $(LibraryFlags) --std=c++20

# By setting this command, all the targets will be executed consequently
install_dependencies: install_portaudio install_fftw
.PHONY: install_dependencies

uninstall_dependencies: uninstall_portaudio uninstall_fftw
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

clean:
	rm -f $(BuildPath)/$(ExecutableName)
.PHONY: clean
