#!/bin/bash

set -x # echo on
set -e # exit on error

brew update
brew install fmt


# Build a specific version of ceres-solver instead of one shipped over brew
curl https://raw.githubusercontent.com/Homebrew/homebrew-core/b0792ccba6e71cd028263ca7621db894afc602d2/Formula/ceres-solver.rb -o ceres-solver.rb
sed -i.bak 's#http://ceres-solver.org/ceres-solver-2.0.0.tar.gz#https://github.com/ceres-solver/ceres-solver/archive/refs/tags/2.1.0rc1.tar.gz#g' ceres-solver.rb
sed -i.bak 's#10298a1d75ca884aa0507d1abb0e0f04800a92871cd400d4c361b56a777a7603#9138a7d80a3142fe3a98519d58a489da9e204b815cd8c0571b3e643f04eca574#g' ceres-solver.rb
sed -i.bak '/patch do/{N;N;N;d;}' ceres-solver.rb
sed -i.bak '/doc.install/d' ceres-solver.rb

brew install --build-from-source ./ceres-solver.rb
