#!/bin/sh
xcodebuild
\cp -r build/Release/libdesignspace.dylib /usr/local/lib/
sudo cp -r build/Release/usr/local/include/. /usr/local/include/designspace/
