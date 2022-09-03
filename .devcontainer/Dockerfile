FROM ros:melodic-perception

# fixed a key issue using the method mentioned in https://github.com/PRBonn/bonnetal/commit/0ab883b8356954c3e57a07a73c30bbac8f035a05
RUN sh -c 'echo "deb http://packages.ros.org/ros/ubuntu $(lsb_release -sc) main" > /etc/apt/sources.list.d/ros-latest.list' && \
  apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-key C1CF6E31E6BADE8868B172B4F42ED6FBAB17C654 && \
  apt update 

ENV CERES_VERSION="1.14.0"
ENV CATKIN_WS=/root/vins_ws

      # set up thread number for building
RUN   if [ "x$(nproc)" = "x1" ] ; then export USE_PROC=1 ; \
      else export USE_PROC=$(($(nproc)/2)) ; fi && \
      apt-get update && apt-get install -y \
      cmake \
      libatlas-base-dev \
      libeigen3-dev \
      libgoogle-glog-dev \
      libsuitesparse-dev \
      python3-catkin-tools \
      git \
      ros-${ROS_DISTRO}-cv-bridge \
      ros-${ROS_DISTRO}-rviz \
      ros-${ROS_DISTRO}-robot-state-publisher \
      ros-${ROS_DISTRO}-image-transport \
      ros-${ROS_DISTRO}-message-filters \
      ros-${ROS_DISTRO}-tf && \
      rm -rf /var/lib/apt/lists/*

      # Build and install Ceres
RUN   git clone https://ceres-solver.googlesource.com/ceres-solver && \
      cd ceres-solver && \
      git checkout tags/${CERES_VERSION} && \
      mkdir build && cd build && \
      cmake .. && \
      make -j$(USE_PROC) install && \
      rm -rf ../../ceres-solver && \
      mkdir -p $CATKIN_WS/src/VINS-Fusion/

# Add VINS-Fusion
WORKDIR $CATKIN_WS/src
RUN git clone https://github.com/HKUST-Aerial-Robotics/VINS-Fusion.git

# Build VINS-Fusion
WORKDIR $CATKIN_WS
RUN /bin/bash -c "source /opt/ros/${ROS_DISTRO}/setup.bash; catkin build; source ${CATKIN_WS}/devel/setup.bash; catkin build"
RUN echo "source ${CATKIN_WS}/devel/setup.bash" >> /root/.bashrc


# change timezone (this is very important otherwise many ROS topic time will be strange)
# TODO: please modify it to your timezone accordingly
ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone


# casadi 

ENV SUPPORT_WS=/root/support_files
WORKDIR $SUPPORT_WS
RUN git clone https://github.com/casadi/casadi.git
RUN cd ${SUPPORT_WS}/casadi && git checkout tags/3.5.5
RUN mkdir ${SUPPORT_WS}/casadi/build && cd ${SUPPORT_WS}/casadi/build
RUN cmake -DWITH_CPLEX=OFF -DWITH_KNITRO=OFF -DWITH_OOQP=OFF -DWITH_SNOPT=OFF ${SUPPORT_WS}/casadi
RUN make -j8
RUN make install

# filter for processing sensor data
WORKDIR $SUPPORT_WS
RUN git clone https://github.com/ShuoYangRobotics/gram_savitzky_golay.git
WORKDIR ${SUPPORT_WS}/gram_savitzky_golay 
RUN git submodule init
RUN git submodule update
RUN mkdir build 
WORKDIR ${SUPPORT_WS}/gram_savitzky_golay/build
RUN cmake -DCMAKE_BUILD_TYPE=Release ../
RUN make 
RUN make install


# for sharing vpg keys
RUN apt-get update && apt-get install gnupg2 -y

# docker build -t vilo_image .

# docker run -d \
# --network host \
# --cap-add=IPC_LOCK --cap-add=sys_nice \
# -v ${PATH_TO_THIS_CODE}/tightly-coupled-visual-inertial-leg-odometry:/root/vins_ws/src/vilo \
# --device /dev/input \
# --name vilo_container \
# vilo_image
