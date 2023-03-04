FROM ros:melodic-perception

ARG USE_PROC=1
RUN echo "USE_PROC=${USE_PROC}"

###################  first system wide configuration ###################################
# change timezone (this is very important otherwise many ROS topic time will be strange)
# TODO: please modify it to your timezone accordingly
ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get update && export DEBIAN_FRONTEND=noninteractive \ 
    && apt-get install -y --no-install-recommends\
    cmake \
    libatlas-base-dev \
    libeigen3-dev \
    libgoogle-glog-dev \
    libsuitesparse-dev \
    python3-catkin-tools \
    lsb-release \
    curl \
    git \
    gdb \
    clang \
    gnupg2 \
    wget \
    openssh-client \
    ros-${ROS_DISTRO}-cv-bridge \
    ros-${ROS_DISTRO}-rviz \
    ros-${ROS_DISTRO}-robot-state-publisher \
    ros-${ROS_DISTRO}-image-transport \
    ros-${ROS_DISTRO}-message-filters \
    ros-${ROS_DISTRO}-tf \
    zsh \
    clang-format-10 && \
    rm -rf /var/lib/apt/lists/*

RUN echo "deb [arch=amd64] http://robotpkg.openrobots.org/packages/debian/pub $(lsb_release -cs) robotpkg" |  tee /etc/apt/sources.list.d/robotpkg.list      
RUN curl http://robotpkg.openrobots.org/packages/debian/robotpkg.key |  apt-key add -

# install python 3.8 but do not break python 2.7 otherwise ros will break
RUN apt-get update && \
    apt-get install -y python3.8 python3.8-dev ipython3
RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.8 1  && \
    curl https://bootstrap.pypa.io/pip/get-pip.py -o get-pip.py  && \
    python3 get-pip.py  && \
    rm get-pip.py  && \
    pip3 --version
###################  configuration specific to the docker user  ###################################
# create a non-root user
ARG USERNAME=EstimationUser
ARG USER_UID=1000
ARG USER_GID=$USER_UID

# Create the user
RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME \
    #
    # [Optional] Add sudo support. Omit if you don't need to install software after connecting.
    && apt-get update \
    && apt-get install -y sudo \
    && echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME \
    && chmod 0440 /etc/sudoers.d/$USERNAME

# [Optional] Set the default user. Omit if you want to keep the default as root.

USER $USERNAME

ENV CERES_VERSION="1.14.0"
ENV CATKIN_WS=/home/${USERNAME}/estimation_ws
ENV SUPPORT_WS=/home/${USERNAME}/support_ws

# Build and install Ceres
RUN mkdir -p $SUPPORT_WS 
WORKDIR $SUPPORT_WS
RUN pwd
RUN git clone https://ceres-solver.googlesource.com/ceres-solver && \
    cd ceres-solver && \
    git checkout tags/${CERES_VERSION} && \
    mkdir build && cd build && \
    cmake .. && \
    make -j ${USE_PROC} && \
    sudo make install

# Add VINS-Fusion
RUN mkdir -p $CATKIN_WS/src
WORKDIR $CATKIN_WS/src
RUN git clone https://github.com/HKUST-Aerial-Robotics/VINS-Fusion.git

# Build VINS-Fusion
WORKDIR $CATKIN_WS
RUN /bin/bash -c "source /opt/ros/${ROS_DISTRO}/setup.bash; catkin build; source ${CATKIN_WS}/devel/setup.bash; catkin build" && \
    echo "source ${CATKIN_WS}/devel/setup.bash" >> /home/${USERNAME}/.bashrc

# casadi 
WORKDIR $SUPPORT_WS
RUN git clone https://github.com/casadi/casadi.git && \
    cd ${SUPPORT_WS}/casadi && git checkout tags/3.5.5 && \
    mkdir ${SUPPORT_WS}/casadi/build && cd ${SUPPORT_WS}/casadi/build && \
    cmake -DWITH_CPLEX=OFF -DWITH_KNITRO=OFF -DWITH_OOQP=OFF -DWITH_SNOPT=OFF ${SUPPORT_WS}/casadi && \
    make -j ${USE_PROC} && \
    sudo make install

# filter for processing sensor data
WORKDIR $SUPPORT_WS
RUN git clone https://github.com/ShuoYangRobotics/gram_savitzky_golay.git && \
    cd gram_savitzky_golay  && \
    git submodule init     && \
    git submodule update   && \
    mkdir build    && \
    cd build   && \
    cmake -DCMAKE_BUILD_TYPE=Release ../   && \
    make -j ${USE_PROC} && \
    sudo make install

# install pinnochio
WORKDIR $SUPPORT_WS
RUN git clone --recursive https://github.com/stack-of-tasks/pinocchio
RUN sudo apt install robotpkg-py27-eigenpy && export CMAKE_PREFIX_PATH=/opt/openrobots
RUN cd pinocchio && mkdir build && cd build && \
    cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local -DBUILD_PYTHON_INTERFACE=OFF && \
    make -j ${USE_PROC}
WORKDIR $SUPPORT_WS/pinocchio/build
RUN sudo make install && \
    echo "export PATH=/usr/local/bin:$PATH"  >> /home/${USERNAME}/.bashrc \ 
    && echo "export PKG_CONFIG_PATH=/usr/local/lib/pkgconfig:$PKG_CONFIG_PATH"  >> /home/${USERNAME}/.bashrc \
    && echo "export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH"  >> /home/${USERNAME}/.bashrc \
    && echo "export PYTHONPATH=/usr/local/lib/python2.7/site-packages:$PYTHONPATH"  >> /home/${USERNAME}/.bashrc \
    && echo "export CMAKE_PREFIX_PATH=/usr/local:$CMAKE_PREFIX_PATH"  >> /home/${USERNAME}/.bashrc

WORKDIR $CATKIN_WS
RUN catkin build
RUN /bin/bash -c "source ${CATKIN_WS}/devel/setup.bash" && echo "source ${CATKIN_WS}/devel/setup.bash" >> /home/${USERNAME}/.bashrc


# add OSQP 
# follow https://osqp.org/docs/get_started/sources.html#build-from-sources to install OSQP from sources
WORKDIR $SUPPORT_WS
RUN git clone --recursive https://github.com/oxfordcontrol/osqp && cd osqp && mkdir build  && cd build && \
    cmake -G "Unix Makefiles" ..  && \
    cmake --build .  && \
    sudo cmake --build . --target install


# add OSQP-python
RUN pip3 install osqp 

# add osqp-eigen
WORKDIR $SUPPORT_WS
RUN git clone https://github.com/robotology/osqp-eigen.git && cd osqp-eigen && mkdir build  && cd build && \
    cmake ../  && \
    make  -j ${USE_PROC} && \
    sudo make install

RUN sh -c "$(wget -O- https://github.com/deluan/zsh-in-docker/releases/download/v1.1.4/zsh-in-docker.sh)" -- \
    -t robbyrussell \
    -a 'source /opt/ros/melodic/setup.zsh' \
    -a 'alias wssetup=". ./devel/setup.zsh"' \
    -p git \
    -p ssh-agent

# config zsh 
CMD ["/bin/zsh"]
ENV SHELL /bin/zsh