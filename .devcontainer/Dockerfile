FROM ros:melodic-perception

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
      python-catkin-tools \
      ros-${ROS_DISTRO}-cv-bridge \
      ros-${ROS_DISTRO}-image-transport \
      ros-${ROS_DISTRO}-message-filters \
      ros-${ROS_DISTRO}-tf && \
      rm -rf /var/lib/apt/lists/* && \
      # Build and install Ceres
      git clone https://ceres-solver.googlesource.com/ceres-solver && \
      cd ceres-solver && \
      git checkout tags/${CERES_VERSION} && \
      mkdir build && cd build && \
      cmake .. && \
      make -j$(USE_PROC) install && \
      rm -rf ../../ceres-solver && \
      mkdir -p $CATKIN_WS/src/VINS-Fusion/

# Add VINS-Fusion
WORKDIR $CATKIN_WS/src
# use the following line if you only have this dockerfile
RUN git clone https://github.com/HKUST-Aerial-Robotics/VINS-Fusion.git

#TODO: add visualizer

# Build VINS-Fusion
WORKDIR $CATKIN_WS
RUN /bin/bash -c "source /opt/ros/${ROS_DISTRO}/setup.bash; catkin build; source ${CATKIN_WS}/devel/setup.bash; catkin build"
RUN echo "source ${CATKIN_WS}/devel/setup.bash" >> /root/.bashrc