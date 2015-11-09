#!/bin/sh

diskutil erasevolume HFS+ 'Ram' `hdiutil attach -nomount ram://8388608`