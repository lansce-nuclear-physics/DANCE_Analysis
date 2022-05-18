#!/bin/bash
kill -9 `ps ux | grep RunAnalysis | awk '{print $2}'`
kill -9 `ps ux | grep DANCE_Analysis | awk '{print $2}'`
