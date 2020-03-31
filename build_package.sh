#!/usr/bin/env bash

R -e 'setwd("/Users/sl31/Documents/barman");devtools::document();devtools::build()'
R -e 'devtools::install()'
R -e 'setwd("/Users/sl31/Documents/barman");pkgdown::build_site();pkgdown::deploy_to_branch(pkg=".")'
