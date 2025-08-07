# Copyright 2025 Giorgio Gamba

# NOTE: $@ stands for the files to build (main.cpp), $^ is for its dependencies

# Configuration variables to be edited as you prefer
BuildPath =  build
ExecutableName = avil

LibraryFlags = -I./lib/portaudio/include ./lib/portaudio/lib/.libs/libportaudio.a -framework CoreAudio -framework AudioToolbox -framework CoreServices -pthread

$(BuildPath)/$(ExecutableName): src/main.cpp
	g++ -o $@ $^ $(LibraryFlags) --std=c++20

install_dependencies:
	mkdir -p lib 
	curl https://files.portaudio.com/archives/pa_stable_v190700_20210406.tgz | tar -zx -C lib
	cd lib/portaudio && ./configure && $(MAKE) -j
.PHONY: install_dependencies

uninstall_dependencies:
	cd lib/portaudio && $(MAKE) uninstall
	rm -rf lib/portaudio
.PHONY: uninstall_dependencies

clean:
	rm -f $(BuildPath)/$(ExecutableName)
.PHONY: clean