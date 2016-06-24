#!/usr/bin/env bash

# it's probably safe to assume for now that if you're using this installer, you're on linux
# for now, all commands assume a python 3.* environment

# python does not have decent XML parsing built in, so this handles object deserialization
pip3 install xmltodict
# this library wraps python's native HTTP library and makes web requests less unpleasant
pip3 install requests
pip3 install dataset # for sqlite

# uncomment and run these instead for systems where pip and python are already assumed as 3.*
#pip install xmltodict
#pip install requests