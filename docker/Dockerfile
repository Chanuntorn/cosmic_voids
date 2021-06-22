FROM ubuntu:20.04

RUN apt -y update && apt install -y apt-utils && echo 2

RUN DEBIAN_FRONTEND=noninteractive \
    apt install -y --no-install-recommends \
    make \
    cmake \
    m4 \
    g++ \
    gcc \
    git \
    file \
    libcairo2 \
    libcairo2-dev \
    libsigc++-2.0-dev \
    libfreetype-dev \
    emacs \
    less \
    patch \
    pkg-config \
    wget \
    libcfitsio-dev \
    libcfitsio-bin \
    libbz2-dev \
    python3 \
    python3-dev \
    python3-pip \
    python3-pil \
    python3-tk \
    # # Remove APT files
    && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Python related stuff
RUN echo "../site-packages" > /usr/local/lib/python3.8/dist-packages/site-packages.pth

# Pip installs
RUN for x in \
    setuptools \
    wheel \
    numpy \
    scipy \
    matplotlib \
    astropy \
    fitsio \
    numba \
    pyyaml \
    requests \
    scikit-learn \
    healpy \
    seaborn \
    photutils \
    ipython \
    ipykernel \
    h5py \
    pytest \
    configobj \
    sqlalchemy \
    corner \
    ; do pip3 install $x; done \
    && rm -Rf /root/.cache/pip

# python = python3
RUN ln -s /usr/bin/python3 /usr/bin/python

RUN mkdir -p /src
WORKDIR /src

RUN git clone https://bitbucket.org/cosmicvoids/vide_public/src/master vide
RUN cd vide \
    && python setup.py build \
    && python setup.py install

RUN git clone https://github.com/Chanuntorn/cosmic_voids cosmic_voids \
    && cd cosmic_voids && git describe --always
RUN cd /src/cosmic_voids && git pull && git describe --always && echo 1

ENV PYTHONPATH /src/cosmic_voids/py:$PYTHONPATH
ENV PATH /src/cosmic_voids/bin:$PATH

RUN mkdir /homedir && chmod 777 /homedir
ENV HOME /homedir

# set prompt and default shell
SHELL ["/bin/bash", "-c"]
ENTRYPOINT ["/bin/bash", "-c"]
CMD ["/bin/bash"]

RUN echo "export PS1='[container] \\u@\\h:\\w$ '" >> $HOME/.bashrc \
  # Create config files in $HOME
  && python -c "import astropy" \
  && python -c "import matplotlib.font_manager as fm; f = fm.FontManager()" \
  && ipython -c "print('hello')" \
  # Download astropy site locations and USNO sky model
  && python -c "from astropy.coordinates import EarthLocation; EarthLocation._get_site_registry(force_download=True)" \
  && python -c "from astropy.coordinates import EarthLocation, SkyCoord, AltAz; from astropy.time import Time; print(EarthLocation.of_site('ctio')); print(SkyCoord(180.,-45.,unit='deg').transform_to(AltAz(obstime=Time(56806.0, format='mjd'), location=EarthLocation.of_site('ctio'))))" \
  # Download astropy IERS leap-second list
  && python -c "from astropy.time import Time; Time.now()" \
  # Make astropy cache files readable!?!!
  && chmod -R a+rwX $HOME/.astropy \
  # Make ipython config files readable!?!!
  && chmod -R a+rwX $HOME/.ipython \
  && echo 1

# Curses upon you, astropy
#RUN wget -O /opt/conda/lib/python3.8/site-packages/astropy/utils/iers/data/Leap_Second.dat \
#    https://raw.githubusercontent.com/astropy/astropy/master/astropy/utils/iers/data/Leap_Second.dat \
#    && rm $HOME/.wget-hsts
