$ipstr=netsh interface ip show address "Wi-fi" | findstr "IP Address"
$iparr = -split $ipstr
$ip=$iparr[-1]
$DISPLAY=":0.0"
$ipfinal=$ip+$DISPLAY

docker run -it --name openfoam4atp  --rm --security-opt seccomp=unconfined -v ~/.Xauthority:/root/.Xauthority -e DISPLAY=$ipfinal -v volOpenFOAM4ATP:/home/openfoam -h openfoam4atp openfoam/openfoam4atp