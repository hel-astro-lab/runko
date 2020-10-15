
for lap in {5400..7400..100}
do
echo "lap is $lap"
    python3 plot3d.py --conf y8x1024s0.ini --lap $lap --var dens
done

#python3 plot3d.py --conf d3x128s10.ini --lap 1000 --var jz
#python3 plot3d.py --conf d3x128s10.ini --lap 1100 --var jz
#python3 plot3d.py --conf d3x128s10.ini --lap 1200 --var jz
#python3 plot3d.py --conf d3x128s10.ini --lap 1300 --var jz
#python3 plot3d.py --conf d3x128s10.ini --lap 1400 --var jz
#python3 plot3d.py --conf d3x128s10.ini --lap 1500 --var jz
#python3 plot3d.py --conf d3x128s10.ini --lap 1600 --var jz
#python3 plot3d.py --conf d3x128s10.ini --lap 1700 --var jz
#python3 plot3d.py --conf d3x128s10.ini --lap 1800 --var jz
#python3 plot3d.py --conf d3x128s10.ini --lap 1900 --var jz
