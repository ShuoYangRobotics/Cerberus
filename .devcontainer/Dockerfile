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


# change timezone (this is very important otherwise many ROS topic time will be strange)
ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone


### run docker as a ssh daemon
RUN apt-get update \
  && apt-get install -y ssh 
RUN touch /etc/ssh/sshd_config_test_clion
RUN ( \
    echo 'LogLevel DEBUG2'; \
    echo 'PermitRootLogin yes'; \
    echo 'PasswordAuthentication yes'; \
    echo 'Subsystem sftp /usr/lib/openssh/sftp-server'; \
    echo 'Port 2237'; \
  ) > /etc/ssh/sshd_config_test_clion \
  && mkdir /run/sshd

# Change the password 'password' to something more secure
RUN useradd -m user && yes password | passwd user
RUN echo 'root:2022quad' | chpasswd

CMD ["/usr/sbin/sshd", "-D", "-e", "-f", "/etc/ssh/sshd_config_test_clion"]


# docker build -t a1_vins_image .

# docker run -d \
# --network host \
# --cap-add=IPC_LOCK --cap-add=sys_nice \
# -v /home/shuoyang/Documents/tightly-coupled-visual-inertial-leg-odometry:/root/vins_ws/src/vilo \
# --device /dev/input \
# --name a1_vins_container \
# a1_vins_image