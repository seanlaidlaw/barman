#!/usr/bin/env bash

R -e 'setwd("/Users/sl31/Documents/barman")'
R -e 'devtools::document();devtools::build()'
R -e 'devtools::install()'
