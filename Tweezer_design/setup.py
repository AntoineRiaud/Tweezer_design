# -*- coding: utf-8 -*-
"""
Created on Sun Apr 03 21:53:14 2016

@author: Antoine Riaud
"""

from distutils.core import setup

setup(name='Acoust_Tweezers_Anisotropic',
      version='2.0',
      description='Acoust_tweezers design package for anisotropic media',
      author='Antoine Riaud, Jean-Louis Thomas, Michael Baudoin, Olivier Bou Matar',
      author_email='antoine.riaud@gmail.com',
      packages=['Tweezer_design'],
      data_files=[('',['Tweezer_design/LiNbO3.mat','Tweezer_design/reticule.svg','Tweezer_design/mirror.svg'])],
     )