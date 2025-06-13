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
ffmpeg -pattern_type glob -i '*.png' -r 100 -b:v 20M -vcodec libopenh264 "AndersonJacksonFoam slowed 10x.mp4"

ffmpeg -pattern_type glob -i '*0.png' -r 100 -b:v 20M -vcodec libopenh264 "AndersonJacksonFoam.mp4"
```