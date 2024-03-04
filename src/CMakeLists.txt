# cmake_minimum_required(VERSION <specify CMake version here>)
project(HgvsGo)

set(CMAKE_CXX_STANDARD 14)

add_executable(HgvsGo main.cpp Variant.cpp Genome.cpp Transcript.cpp Util.cpp Hgvsc.cpp SeqUtil.cpp Hgvsp.cpp Translator.cpp Mrna.cpp Hgvs.cpp TranscriptIntervalForest.cpp)