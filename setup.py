#!/usr/bin/env python

from distutils.core import setup


jbastro_dep_link='https://github.com/baileyji/jbastro/tarball/v0.3#egg=jbastro-0.3'

setup(name='M2FS',
      version='0.4',
      description='M2FS reduction utilities',
      author='Jeb Bailey',
      author_email='baileyji@umich.edu',
      url='',
      packages=['m2fs','m2fs.obs','m2fs.plate'],
      scripts=['bin/m2fs_merge.py',
               'bin/m2fs_obslog.py',
               'bin/m2fs_platesummary.py',
               'bin/m2fs_scatter.py',
               'bin/m2fs_stack.py',
               'bin/m2fs_focus.py'],
      install_requires=['numpy','jbastro>=0.3'],
      dependency_links = [jbastro_dep_link]
      )
