cmake_minimum_required(VERSION 3.10)
project(impl C)

# Znajdź wszystkie pliki .c w bieżącym folderze
file(GLOB SOURCES "*.c")

# Stwórz executable z wszystkich znalezionych plików
add_executable(program ${SOURCES})