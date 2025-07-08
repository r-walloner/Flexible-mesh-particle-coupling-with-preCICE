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
ffmpeg -pattern_type glob -framerate 500 -i "*.png" -vcodec libopenh264 -b:v 25M -r 30 "../$(basename "$(pwd)").mp4" && \
ffmpeg -pattern_type glob -framerate 250 -i "*.png" -vcodec libopenh264 -b:v 25M -r 30 "../$(basename "$(pwd)"), slowed 5x.mp4" && \
ffmpeg -pattern_type glob -framerate 50 -i "*.png" -vcodec libopenh264 -b:v 25M -r 30 "../$(basename "$(pwd)"), slowed 10x.mp4"
```


Generate `0/alphaSolid`
1. Generate run.
2. In `system/controlDict`, set `writeInterval` to equal `deltaT` and disable `writeCompression`.
3. Run for two timewindows.
4. Copy `0.00002/alphaSolid` to `0/alphaSolid`.
