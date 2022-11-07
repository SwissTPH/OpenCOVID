#!/usr/bin/env python
# -*- coding: utf-8 -*-

import yaml

# ---------------------------------------------------------
# Read yaml file using 'yaml' from 'pyymal' package
# ---------------------------------------------------------
def read_yaml_py(file):
  
  # Open connection and parse yaml syntax
  with open(file, 'r') as stream:
      y = yaml.safe_load(stream)
  
  return y
