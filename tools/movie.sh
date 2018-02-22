#general with 0001, 0002, etc naming convention
ffmpeg -f image2 -r 40 -i movie/filename_%03d.png -c:v libx264 -pix_fmt yuv420p movie/name.mp4

#or using pattern_type that is more flexible
#ffmpeg -f image2 -r 15 -pattern_type glob -i 'out/vlasov*.png' -c:v libx264 -pix_fmt yuv420p out/twostream.mp4
ffmpeg -f image2 -r 15 -pattern_type glob -i 'node_*.png' -c:v libx264 -pix_fmt yuv420p node.mp4
