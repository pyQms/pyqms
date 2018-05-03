#!/bin/bash

echo 'Building package...'
if [ ! -d "dist" ]; then
  mkdir dist
  echo 'Created directory dist'
else
  rm -rf dist/*
  echo 'Removed previous distribution from dist'
fi

rm -rf docs/build/*

# Evoke Sphinx to create html and pdf documentation
cd docs
make html latexpdf
cd ..


# Creating Python packages
python3 setup.py sdist --formats=bztar,gztar,zip
# python3.4 setup.py sdist --formats=zip
cd dist
tar xvfj *.bz2
cd ..
