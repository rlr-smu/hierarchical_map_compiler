#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <iostream>
#include <string>

extern std::string project_script_path;

int main(int argc, char *argv[]) {
  if (argc < 2) {
    std::cerr << "Usage: <path_to_script_directory>" << std::endl;
    exit(1);
  }
  project_script_path = argv[1];
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
