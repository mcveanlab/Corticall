#!/bin/bash

set -euxo pipefail

cd wdl/; rm -f wdls.zip; zip -r wdls.zip * > /dev/null; cd ..;

#cromshell submit wdl/ProcessPfCross.wdl input/ProcessPfCross/ProcessPfCross.3D7xHB3.json options/default.json wdl/wdls.zip
#cromshell submit wdl/ProcessPfCross.wdl input/ProcessPfCross/ProcessPfCross.HB3xDD2.json options/default.json wdl/wdls.zip
#cromshell submit wdl/ProcessPfCross.wdl input/ProcessPfCross/ProcessPfCross.7G8xGB4.json options/default.json wdl/wdls.zip
cromshell submit wdl/ProcessPfCross.wdl input/ProcessPfCross/ProcessPfCross.803xGB4.json options/default.json wdl/wdls.zip

#cromshell submit wdl/ProcessPfCross.wdl input/ProcessPfCross/ProcessPfCross.test.json options/default.json wdl/wdls.zip
