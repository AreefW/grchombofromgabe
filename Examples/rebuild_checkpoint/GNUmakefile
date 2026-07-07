# -*- Mode: Makefile -*-
#
# Minimal Chombo GNUmakefile for rebox_checkpoint_level0.cpp.
#
# Put this file in a small exec directory with rebox_checkpoint_level0.cpp,
# rename it to GNUmakefile, and set CHOMBO_HOME to your Chombo lib directory.
makefiles += rebox_checkpoint_level0
# Executable basename. Chombo's Make.example expects rebox_checkpoint_level0.cpp
# in one of src_dirs and builds an executable with this basename plus the usual
# Chombo dimension/MPI/debug suffixes.
ebase = rebox_checkpoint_level0
# Adjust this path. It must point at the Chombo "lib" directory, not the repo root.
# CHOMBO_HOME = ../../Chombo/lib
# AMRIO, LevelData, domainSplit, DisjointBoxLayout, and LoadBalance live in the
# BoxTools layer in the GRTLCollaboration/Chombo tree. Keeping AMRTools here is
# harmless in normal Chombo builds and helps with installations that put AMR
# hierarchy helpers there.
LibNames = AMRTools BoxTools
base_dir = .
src_dirs = 
include $(CHOMBO_HOME)/mk/Make.example