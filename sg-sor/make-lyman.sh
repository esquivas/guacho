#!/bin/bash

echo C0
cp lyman_alpha_tau-c0.f90 lyman_alpha_tau.f90
make lyman_alpha_tau
mv lyman_alpha_tau lyman_alpha_tau-c0

echo H0
cp lyman_alpha_tau-h0.f90 lyman_alpha_tau.f90
make lyman_alpha_tau
mv lyman_alpha_tau lyman_alpha_tau-h0

echo C02T
cp lyman_alpha_tau-c02T.f90 lyman_alpha_tau.f90
make lyman_alpha_tau
mv lyman_alpha_tau lyman_alpha_tau-c02t

echo H0
cp lyman_alpha_tau-h02T.f90 lyman_alpha_tau.f90
make lyman_alpha_tau
mv lyman_alpha_tau lyman_alpha_tau-h02t
