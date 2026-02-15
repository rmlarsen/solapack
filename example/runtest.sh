#! /bin/csh -f
setenv PLAT `grep PLAT ../make.inc | cut -d" " -f 3`
echo  ********* Setting up mode kernels **********
../get-config set-2drls.cfg | ../set-2drls.nn.${PLAT}.x
echo  ********* Performing inversion **********
../get-config 2dsola.cfg | ../2dsola.lanczos.${PLAT}.x

