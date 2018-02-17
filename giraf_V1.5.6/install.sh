#!/bin/ksh
#programme installation GIRAF
# E. Chaxel - Atmo Rhone-Alpes - 2010

localdir=`pwd`

if [ "${F90}" == "" ] ; then
  echo "Indiquer le compilateur FORTRAN 90 en declarant"
  echo "export F90=nom_du_compilateur (ex: gfortran, ifort)"
  exit
else
  echo "Utilise F90="${F90}
fi

if [ "${NETCDFHOME}" == "" ] ; then
  echo "Indiquer l emplacement de la bibliotheque NetCDF en declarant"
  echo "export NETCDFHOME=netcdf_root (ex: /opt/netcdf)"
  exit
else
  echo "Utilise NETCDFHOME="${NETCDFHOME} 
fi

if [ "${PROJECTIONLIB}" == "" ] ; then
  echo "Indiquer l emplacement de la bibliotheque PROJECTION en declarant"
  echo "export PROJECTIONLIB=projection (ex: /opt/projection)"
  exit
else
  echo "Utilise PROJECTIONLIB="${PROJECTIONLIB} 
fi

# nom de l executable (ne pas toucher)
exe=giraf.exe

# REGLE SPECIALE GFORTRAN
case ${F90} in 
gfortran)
export FCOPT="-O2 -ffree-form -ffree-line-length-none -fbounds-check";;
*)
export FCOPT="-O2";;
esac

# options de compilation personnalisée (DECOMMENTER LES LIGNES SUIVANTES)
#export F90="gfortran"
#export FCOPT="-O2 -ffree-form -ffree-line-length-none -fbounds-check"
#export NETCDFHOME="/opt/netcdf-4.0.1-gfortran-64"

# NE PLUS MODIFIER LA SUITE...

case $1 in 
clean)
rm -rf ${exe}
cd src
make clean
cd ..
echo "Desinstallation du mailleur complete. Au revoir."
;;
*)
cd src
echo "Programme en cours de compilation..."
make clean >  ${localdir}/make.log 2>&1
make       >> ${localdir}/make.log 2>&1
cd ..
if [ -f src/${exe} ] ; then
  ln -sf src/${exe} .
  echo "Installation du raffinement de maillage complete. Merci."
  echo "Pour de l aide, utiliser ${exe} -h"
else
  echo "Erreur dans la compilation: jeter un oeil au fichier make.log"
  exit
fi
;;
esac
