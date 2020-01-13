C:\ffmpeg\bin\ffmpeg.exe -r 5 -f image2 -i "soc_%%06d.png" -crf 1 -q:v 1 -s 1000x800 -vcodec wmv2 _vid.wmv
PAUSE