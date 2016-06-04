#pragma once

// Ensight file output
// case file
// geometry file - domain, then faces in groups by type/connectivity
// per-element, per-node, etc. files.


struct EnsightOut
{
	EnsightOut()
	{
	
	}
};

/*
FORMAT
type:  ensight gold
GEOMETRY
model:  Splitterplate.geo
VARIABLE
scalar per node:	 r	r
scalar per node:	 P	P
vector per node:	 v	v
scalar per node:	 u	u
scalar per node:	 m	m
scalar per node:	 t	t
scalar per node:	 Terror	Terror
scalar per node:	 Perror	Perror
scalar per element:	 yplus	yplus
scalar per element:	 matchyplus matchyplus
*/