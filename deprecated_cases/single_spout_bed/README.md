LIGGGHTS CMake options:
```
ENABLE_COHESION_OFF
ENABLE_HOOKE_STIFFNESS
ENABLE_ROLLING_OFF
ENABLE_SURFACE_DEFAULT
ENABLE_TANGENTIAL_HISTORY
```

Render pngs with ffmpeg:
```shell
ffmpeg -pattern_type glob -framerate 100 -i "*.png" -vcodec libopenh264 -b:v 25M -r 30 "../../$(basename "$(pwd)"), slowed 10x.mp4"

ffmpeg -pattern_type glob -framerate 1000 -i "*.png" -vcodec libopenh264 -b:v 25M -r 30 "../../$(basename "$(pwd)").mp4"
```