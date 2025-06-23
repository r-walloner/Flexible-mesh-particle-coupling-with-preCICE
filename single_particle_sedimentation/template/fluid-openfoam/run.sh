#!/usr/bin/env bash
set -e -u

. ../../../../other/tools/log.sh
exec > >(tee --append "$LOGFILE") 2>&1

blockMesh

../../../../other/tools/run-openfoam.sh "$@"
. ../../../../other/tools/openfoam-remove-empty-dirs.sh && openfoam_remove_empty_dirs

close_log
