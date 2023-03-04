XAUTH=/tmp/.docker.xauth
touch $XAUTH
if [ ! -f $XAUTH ]
then
    xauth_list=$(xauth nlist :0 | sed -e 's/^..../ffff/')
    if [ ! -z "$xauth_list" ]
    then
        sudo echo $xauth_list | sudo xauth -f $XAUTH nmerge -
    else
        sudo touch $XAUTH
    fi
    sudo chmod a+r $XAUTH
fi

