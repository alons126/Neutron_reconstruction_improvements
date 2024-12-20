#!/bin/bash

# Path to the directory containing the source files
SOURCE_DIR="/Users/alon/Projects/Neutron_reconstruction_improvements"
echo SOURCE_DIR
echo 

# Output CMakeLists.txt file
CMAKE_FILE="${SOURCE_DIR}/CMakeLists.txt"
echo CMAKE_FILE
echo 

rm -rf CMAKE_FILE

# Create or overwrite the CMakeLists.txt file
echo "cmake_minimum_required(VERSION 3.16)" > "$CMAKE_FILE"
echo "project(Neutron_reconstruction_improvements LANGUAGES CXX)" >> "$CMAKE_FILE"
echo "" >> "$CMAKE_FILE"
echo "set(CMAKE_CXX_STANDARD 17)" >> "$CMAKE_FILE"
echo "set(CMAKE_CXX_STANDARD_REQUIRED ON)" >> "$CMAKE_FILE"
echo "" >> "$CMAKE_FILE"

# Find all `.cpp` files in the source directory
echo "set(SOURCES" >> "$CMAKE_FILE"
find "$SOURCE_DIR" -name "*.cpp" -print0 | while IFS= read -r -d '' file; do
    # Append each file, properly quoted, to the SOURCES list
    echo "    \"$file\"" >> "$CMAKE_FILE"
done
echo ")" >> "$CMAKE_FILE"
echo "" >> "$CMAKE_FILE"

# Add the executable definition
echo "add_executable(\${PROJECT_NAME} \${SOURCES})" >> "$CMAKE_FILE"

echo "CMakeLists.txt has been updated!"


# echo "set(SOURCES" > "$output_file"
# find "$src_dir" -type f \( -name "*.cpp" -o -name "*.c" -o -name "*.C" \) | while read -r file; do
#     echo "    $file" >> "$output_file"
# done
# echo ")" >> "$output_file"

# echo "add_executable(Neutron_reconstruction_improvements \${SOURCES})" >> "$output_file"
