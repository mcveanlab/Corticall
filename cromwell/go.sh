#!/bin/bash

set -euxo pipefail

cd wdl/; rm -f wdls.zip; zip -r wdls.zip * > /dev/null; cd ..;

#cromshell submit wdl/AnnotateVCFs.wdl input/AnnotateVCFs/AnnotateVCFs.3D7xHB3.json options/default.json wdl/wdls.zip
#cromshell submit wdl/AnnotateVCFs.wdl input/AnnotateVCFs/AnnotateVCFs.HB3xDD2.json options/default.json wdl/wdls.zip
#cromshell submit wdl/AnnotateVCFs.wdl input/AnnotateVCFs/AnnotateVCFs.7G8xGB4.json options/default.json wdl/wdls.zip
#cromshell submit wdl/AnnotateVCFs.wdl input/AnnotateVCFs/AnnotateVCFs.803xGB4.json options/default.json wdl/wdls.zip

cromshell submit wdl/Simulate.wdl input/Simulate/Simulate.json options/default.json wdl/wdls.zip
