#!/bin/bash
# Helpful Homebrew script for installing commonly-used Mac OSX apps.

ruby -e "$(curl -fsSL https://raw.github.com/Homebrew/homebrew/go/install)" # download brew
brew doctor

brew install git
brew install homebrew/versions/gcc49
brew install python3
brew install aspell
