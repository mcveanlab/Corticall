#!/bin/bash

set -euxo pipefail

cd wdl/; rm -f wdls.zip; zip -r wdls.zip * > /dev/null; cd ..;

cromshell submit wdl/ProcessPfCross.wdl input/ProcessPfCross/ProcessPfCross.test.json options/default.json wdl/wdls.zip

#cromshell submit wdl/PfAlign.wdl input/PfAlign/PfAlign.all_crosses.json options/default.json wdl/wdls.zip

#cromshell submit wdl/PfAsm.wdl input/PfAsm/PfAsm.test.json options/default.json wdl/wdls.zip
#cromshell submit wdl/PfAsm.wdl input/PfAsm/PfAsm.3D7xHB3.json options/default.json wdl/wdls.zip
#cromshell submit wdl/PfAsm.wdl input/PfAsm/PfAsm.HB3xDD2.json options/default.json wdl/wdls.zip
#cromshell submit wdl/PfAsm.wdl input/PfAsm/PfAsm.7G8xGB4.json options/default.json wdl/wdls.zip
#cromshell submit wdl/PfAsm.wdl input/PfAsm/PfAsm.803xGB4.json options/default.json wdl/wdls.zip
