#!/bin/sh
set -e

mkdir -p "${PREFIX}/bin"
mkdir -p "${PREFIX}/docs"
mkdir -p "${PREFIX}/accessory"
mkdir -p "${PREFIX}/bin/modules"

cp -r bin/* "${PREFIX}/bin/"
cp -r bin/modules/* "${PREFIX}/bin/modules"
cp -r docs/* "${PREFIX}/docs/"
cp -r accessory/* "${PREFIX}/accessory/"

