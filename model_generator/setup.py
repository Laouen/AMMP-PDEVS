#!/usr/bin/env python

from setuptools import setup, find_packages

setup(name='PMGBP',
      version='0.1',
      description='Parsing and Model Generation of Biological Processes',
      author='Laouen M.L. Belloli',
      author_email='laouen.belloli@gmail.com',
      packages=find_packages(),
      scripts=['scripts/pmgbp_generate_model'],
      package_data={'PMGBP': ['templates/atomic_model_definition.tpl.hpp', 'templates/dynamic_atomic.tpl.hpp', 'templates/dynamic_coupled.tpl.hpp', 'templates/dynamic_defined_atomic.tpl.hpp', 'templates/dynamic_reaction_set.tpl.hpp']},
      install_requires=['python-gflags', 'BeautifulSoup', 'lxml', 'tqdm'])