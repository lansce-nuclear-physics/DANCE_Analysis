#!/bin/bash

kill -9 `ps aux | grep DANCE_Analysis | awk '{print $2}'`