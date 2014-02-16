#!/usr/bin/env python

from distutils.core import setup

setup(name='M2FS',
      version='0.2',
      description='M2FS reduction utilities',
      author='Jeb Bailey',
      author_email='baileyji@umich.edu',
      url='',
      packages=['m2fs','m2fs.obs','m2fs.plate'],
      scripts=['bin/m2fs_merge.py',
               'bin/m2fs_obslog.py',
               'bin/m2fs_platesummary.py',
               'bin/m2fs_scatter.py',
               'bin/m2fs_stack.py'],
      install_requires=['numpy']
      )