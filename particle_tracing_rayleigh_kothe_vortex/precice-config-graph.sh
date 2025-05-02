#!/bin/bash

precice-config-visualizer \
 --data-access full ./precice-config.xml \
 --communicators hide \
 --cplschemes hide \
 --mappings full \
 > precice-config.dot \
&& dot -Tpng -Gdpi=300 -O precice-config.png precice-config.dot
