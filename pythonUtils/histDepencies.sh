#!/bin/bash
echo Downloading pip
sudo apt-get update && sudo apt-get -y upgrade
sudo apt-get install python-pip
echo upgrading pip
pip install --upgrade pip
echo Downloading numpy
pip install --user numpy
echo Downloading matplotlib
pip install --user matplotlib


